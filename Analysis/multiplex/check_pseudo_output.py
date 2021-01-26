import pandas as pd
ad_p = pd.read_csv("../PBMC_P_cellSNP/cellSNP.tag.AD.mtx",sep="\t", comment="%", header=None)
ad_all = pd.read_csv("pseudo/numC1000_ispropFalse/cellSNP.tag.AD.mtx",sep="\t", comment="%", header=None)