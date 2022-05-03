import numpy as np
from glmnet import ElasticNet
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

def fit_clf(gene_gt, expr, clf, tr_indiv=None, indiv_age_map=None, 
            sep=False, sex=None):
    output = {}
    X = np.transpose(gene_gt[tr_indiv].to_numpy())
    coef = []
    perf = [] # Contains r2 performance of age and genetic factors
    if indiv_age_map != None:
        ageX = np.expand_dims(
                             np.array([indiv_age_map[ind] for ind in tr_indiv]), 
                             axis=1)
        if sex: # Add age x sex interaction term as last predictor
            sexX = ageX*np.expand_dims(covar.loc["sex", tr_indiv], axis=1)
            X = np.concatenate((X, sexX), 1)
        if not sep: # Add age as last predictor if age_reg
            X = np.concatenate((X, ageX),1)
        else:
            ageclf = LinearRegression()
            ageclf.fit(ageX, np.transpose(expr[tr_indiv]))
            coef += [ageclf.coef_[0]]
        
    clf.fit(X, np.transpose(expr[tr_indiv]))
    coef = list(clf.coef_) + coef
    return clf.cv_mean_score_[clf.lambda_max_inx_], coef

# Given trained model weights, evaluates R^2 age and genetics
def eval_modelperf(gene_gt, coefs, expr_row, t_indiv, indiv_age_map, sex):
    centered_expr = expr_row[t_indiv]-np.mean(expr_row[t_indiv])
    last_index = len(coefs) if indiv_age_map==None else (-1 if not sex else -2)
    pred = np.matmul(coefs[:last_index], gene_gt[t_indiv].to_numpy())
    pred = pred - np.mean(pred)

    # Variance explained or R^2 for genotype weights in lin model 
    # (after correcting for covariates)
    gene_r2 = r2_score(centered_expr, pred)  

    age_r2 = sex_r2 = None
    # Variance explained by age weights
    if indiv_age_map != None:
        age_r2 = r2_score(
                centered_expr, 
                np.array([indiv_age_map[ind] for ind in t_indiv])*coefs[-1]) 

    # Variance explained by age x sex interaction term
    if indiv_age_map != None and sex:
        sex_r2 = r2_score(
                centered_expr, 
                np.array([indiv_age_map[ind] for ind in t_indiv])\
                        *covar.loc["sex", t_indiv]*coefs[-2])  

    return gene_r2, age_r2, sex_r2


# Run train age/genetic model using bootstrap samples of individuals
def bootstrap_train(gene_gt, row, indiv, tr_indiv, t_indiv=None, 
                    indiv_age_map=None, sep=False, bootstrap=None,
                    sex=False, seed=Seed(0,0), old=True):
    bootres = {}
    prefix = "" if old else "young_"

    for i in range(bootstrap): # repeat analysis n times and sample w replacement
        clf = ElasticNet(alpha=0.5, n_splits=10, n_jobs=1, fit_intercept=True, 
                random_state=seed.get(row["gene_id"], i, t="permute", old=old)[0])

        if bootstrap!=1:
            np.random.seed(seed.get(t="shuffle", old=old))
            tr_indivboot = np.random.choice(tr_indiv, len(tr_indiv))
            t_indivboot = t_indiv
            if tr_indiv==t_indiv: 
                t_indivboot=tr_indivboot
        else:
            tr_indivboot = tr_indiv
            t_indivboot = t_indiv

        t_r2, coef = \
            fit_clf(gene_gt, row, clf, tr_indivboot, 
                    indiv_age_map=indiv_age_map, sep=sep, sex=sex)

        # evaluate r2 of gt and age
        gene_r2, age_r2, sex_r2 = eval_modelperf(gene_gt, coef, row, 
                                         t_indivboot, indiv_age_map, sex=sex)

    bootres.setdefault(prefix+"t_r2",[]).append(t_r2)
    bootres.setdefault(prefix+"gene_r2",[]).append(gene_r2)
    if indiv_age_map != None:
        bootres.setdefault(prefix+"age_r2",[]).append(age_r2)
        bootres.setdefault(prefix+"age_coef",[]).append(coef[-1])
        if sex:
            bootres.setdefault(prefix+"sex_coef",[]).append(coef[-2])
    return bootres

def train_model(gt, row, covar, indiv, tr_old_indiv, tr_young_indiv=None, 
                t_old_indiv=None, t_young_indiv=None, eval_other=False, 
                indiv_age_map=None, reg_out=True, sep=False, bootstrap=None,
                intermresdir=None, sex=False, seed=Seed(0,0)):
    output = {}
    gene_gt = gt.loc[(gt["Chrom"]==row["#chr"]) & 
                     (gt["Pos"]>=(row["start"]-5e5)) & 
                     (gt["Pos"]<=(row["start"]+5e5)), 
                     gt.columns.difference(["Pos", "Chrom"])] # indiv x positions

    if gene_gt.shape[0] == 0: # No variants to be tested, terminate early and results should be 0's
        print("Zero variants for gene %s" % row["gene_id"])
        output["invalid_gene"] = [row["gene_id"]]
        return output # return zero performance
    else:
        output["gene"] = [row["gene_id"]]

    if reg_out: # Regress the covariates out from expression data
        print(indiv)
        skclf = LinearRegression()
        covar_T = np.transpose(covar[indiv].to_numpy())
        expr_T = np.transpose(row[indiv].to_numpy())
        skclf.fit(covar_T, expr_T)
        covar_r2 = skclf.score(covar_T, expr_T)
        row[indiv] = row[indiv] - skclf.predict(covar_T)
        output["covar_r2"] = [covar_r2]

    bootres = bootstrap_train(gene_gt, row, indiv, tr_old_indiv, 
                              t_indiv=t_old_indiv, 
                              indiv_age_map=indiv_age_map, sep=sep,
                              bootstrap=bootstrap, sex=sex, seed=seed, 
                              old=True)
    #output["coef"] = [old_coef]
    for k, v in bootres.items():
        output[k] = [np.mean(bootres[k])]
        output["sd_"+k] = [np.std(bootres[k])]
    if indiv_age_map!=None:
        output["nonz_agecoef"] = [np.mean(bootres["age_coef"]==0)]

    # Build model and evaluate on young cohort if exists
    if tr_young_indiv!=None: # Evaluates r2 on other cohort
        bootres = bootstrap_train(gene_gt, row, indiv, tr_young_indiv, 
                                  t_indiv=t_young_indiv, 
                                  indiv_age_map=indiv_age_map, sep=sep,
                                  bootstrap=bootstrap, sex=sex, seed=seed, 
                                  old=False)
        #output["coef"] = [old_coef]
        for k, v in bootres.items():
            output[k] = [np.mean(bootres[k])]
            output["sd_"+k] = [np.std(bootres[k])]
        if indiv_age_map!=None:
            output["nonz_young_agecoef"] = [np.mean(bootres["young_age_coef"]==0)]

    if intermresdir != None:
        np.save(output, intermresdir+"/"+gt["gene"])
    return output
