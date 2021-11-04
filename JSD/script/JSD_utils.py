import os
import pandas as pd
from scipy.stats import entropy, ttest_ind_from_stats
import numpy as np
from scipy.spatial import distance as dist
import time

def gmt_line_splitter(line:str):
    lst = line.split("\t")
    return lst[0], lst[1], lst[3:]

def gene_mapping(gmt, df:pd.DataFrame):
    "Retuen dictionary of geneset name mapped to corresponding rows from df"
    sets = {}
    with open(gmt, 'r', encoding='utf-8') as genesets:
        for line in genesets.readlines():
            tu = gmt_line_splitter(line)
            sets[tu[0]] = (set(tu[2]), tu[1])
    retdic = {}
    for se in sets.keys():
        part = df.loc[contained_geneset(df, sets[se][0])]
        retdic[se] = (part, sets[se][1], sets[se][0])
    return retdic

def gene_mapping_table(table, df):
    dic = {}
    somelst = [str(i) for i in list(range(1, 23))]
    somelst.extend(['Y', 'X', 'autosome'])
    for sex in somelst[::-1]:
        if sex == 'autosome':
            d = table.loc[table['chromosome_name'] != 'X']['Name.Description']
            e = [val.split('#')[1] for val in d]
            dic[sex] = df.loc[[val.split('#')[1] in e for val in df["Name#Description"].values]]
            continue
        #print(df["Name#Description"].values)
        #print(table.loc[table['chromosome_name'] == sex]['Name.Description'])
        c = [val.split('#')[1] for val in  table.loc[table['chromosome_name'] == sex]['Name.Description']]
        b = [val.split('#')[1] in c for val in df["Name#Description"].values]
        #print(b)
        a = df.loc[b]
        #print(a.values)
        if len(a.index) == 0:
            continue
        dic[sex] = a
    return dic

def pairwise_JSD(df):
    "Return the average of all pairwise JSD saves values of pairwise JSD values to csv file in specified path"
    #print("df size", df.shape)
    df = df.T
    X = df.values
    X = np.transpose(df.values)
    X = (X.T - np.amin(X, 1)).T
    X = (X.T / np.sum(X, 1)).T
    X[np.isnan(X)] = 0
    X = X.T
    d = lambda u, v: entropy((u + v) / 2.0) - 0.5 * (entropy(u) + entropy(v))
    M = dist.pdist(X, d)
    M = np.sqrt(M)
    M[np.isnan(M)] = 0
    return np.mean(M), M

def pairwise_JSD_mask(df, mask):
    "Return the average of masked pairwise JSD saves values of pairwise JSD values to csv file in specified path"
    #print("df size", df.shape)
    print(df)
    mask = mask[np.triu_indices(mask.shape[0], 1)]
    X = df.values
    X = np.transpose(df.values)
    X = (X.T - np.amin(X, 1)).T
    X = (X.T / np.sum(X, 1)).T
    X[np.isnan(X)] = 0
    d = lambda u, v: entropy((u + v) / 2.0) - 0.5 * (entropy(u) + entropy(v))
    M = dist.pdist(X, d)
    M = np.multiply(M, mask)
    M = np.sqrt(M)
    M[np.isnan(M)] = 0
    return np.mean(M), M

def age_search(id_row, table):
    """Returns age row based on the ids"""
    table = table.set_index("SUBJID")
    lst = (table.loc[id_row, "AGE"])
    return np.array(lst)

def split_age(df, table, threshold):
    """Return split of temp dataframe given table of sample info according to threshold"""
    agelist = age_search(df.columns.values[1:], table)
    row = None
    df = df.drop(df.columns[0], axis=1)
    df = df.transpose()
    young = df.loc[agelist < threshold]
    old = df.loc[agelist >= threshold]
    return young, old, row

def split_age_new(df, table, threshold):
    """Return split of temp dataframe given table of sample info according to threshold"""
    young = df[set(table.loc[table['AGE'] < threshold]['SAMPID']).intersection(set(df.columns[1:]))]
    old = df[set(table.loc[table['AGE'] >= threshold]['SAMPID']).intersection(set(df.columns[1:]))]
    return young, old

def split_age_bins(df, table, num_bins):
    """Return binned dataframe"""
    i = 20
    upper = 80
    width = upper // num_bins
    binned_frame = []
    while (i < upper):
        more = table['AGE'] >= i
        less = table['AGE'] < (i + width)
        binned_frame.append(df[set(table.loc[more & less]['SAMPID']).intersection(set(df.columns[1:]))])
        print(df[set(table.loc[more & less]['SAMPID']).intersection(set(df.columns[1:]))])
        i += width
    return binned_frame

def pairwise_JSD_age(df, table):
    "Return the average of all pairwise JSD, pairwise average age and its distribution "
    #print("df size", df.shape)
    print(df)
    table = table.set_index("SAMPID")
    agelist = table.loc[df.columns, "AGE"]
    agelist = np.array([[i] for i in agelist])
    # print(agelist.shape)
    # print(agelist)
    ave = lambda u,v: (u + v) / 2
    pair_age = dist.pdist(agelist, ave)
    print(pair_age)

    df = df.T
    X = df.values
    X = np.transpose(df.values)
    X = (X.T - np.amin(X, 1)).T
    X = (X.T / np.sum(X, 1)).T
    X[np.isnan(X)] = 0
    X = X.T
    d = lambda u, v: entropy((u + v) / 2.0) - 0.5 * (entropy(u) + entropy(v))
    M = dist.pdist(X, d)
    M = np.sqrt(M)
    M[np.isnan(M)] = 0
    print(M, pair_age)
    return np.mean(M), M, pair_age
