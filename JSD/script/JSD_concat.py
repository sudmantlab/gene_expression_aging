import os
import pandas as pd
import glob

# def concat(in, out):
#     path = out[:out.rfind("/")]
#     data = None
#     for root, dirs, files in os.walk(path):
#         for file in files:
#             if data is None:
#                 data = pd.read_csv(root + '/' + file, index_col=0)
#             elif "JSD_result.csv" in file:
#                 frame = pd.read_csv(root + '/' + file, index_col=0)
#                 data = pd.concat([data, frame], ignore_index=True)
#     data.to_csv(out)

def concat(path, extension):
    os.chdir(path)
    files = [file for file in glob.glob('*{}'.format(extension))]
    combined_csv = pd.concat([pd.read_csv(f) for f in files])
    combined_csv.to_csv("All_{}".format(extension), index=False)

def concat_dist(path, extension):
    os.chdir(path)
    files = [file for file in glob.glob('*{}'.format(extension))]
    combined_csv = pd.concat([pd.read_csv(f, header=None).rename({0:f[:f.rfind('_d')]}\
        , axis=1) for f in files], axis=1)
    combined_csv.to_csv("All_{}".format(extension), index=False)

def main():
    # path = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/JSD_All/"
    # extension = "_JSD_result.csv"

    # concat(path, extension)

    path = "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/JSD_geneset/JSD_All/"
    extension = "_dist.csv"
    concat_dist(path, extension)

main()
