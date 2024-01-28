#' @title py_qPCR_analysis
#' @author GM.W
#' @date 2023-05-19
import argparse
import pandas as pd
import numpy as np
from plotnine import (ggplot, aes, geom_jitter, geom_bar, geom_errorbar,
                      theme_bw, facet_wrap, theme, labs, ggsave,
                      element_blank, element_rect, element_text,
                      scale_fill_manual, scale_color_manual)

# parameters
parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', type=str, help='The name of input qPCR file.',required=True)
parser.add_argument('--output', '-o', type=str, help='The name of output file.', required=True)
parser.add_argument('--control', '-c', type=str, help='The name of control in experiment.', required=True)
parser.add_argument('--reference', '-r', type=str, help='The reference gene name.', required=True)
parser.add_argument('--fold', '-f', type=float, help='The foldchange of sd to deal with outliar values. Default as 3.', default=3)

args = vars(parser.parse_args())

infile = args["input"]
outfile = args["output"]
ref = args["reference"]
con = args["control"]
fold = args["fold"]

# read data
print("Reading data...")
informat = infile[-3:]
match informat:
    case "csv":
        data = pd.read_csv(infile)
    case "xls":
        data = pd.read_excel(infile)
    case "lsx":
        data = pd.read_excel(infile)
    case "txt":
        data = pd.read_table(infile, sep="\t")
    case _:
        print("Not supported format.")

print("Analyzing...")
data.columns = ["samples","gene","ct"]

sample_list = data["samples"].drop_duplicates().values
gene_list = data["gene"].drop_duplicates().values

# mean & sd
for i in sample_list:
    for j in gene_list:
        select = (data["samples"]==i) & (data["gene"] == j)
        data.loc[select,"ct_mean"] = data.loc[select,"ct"].mean()
        data.loc[select,"ct_sd"] = data.loc[select,"ct"].std()

# filter outliars
for i in range(data.shape[0]):
    data.loc[i,"keep"] = False if np.abs(data.loc[i,"ct"] - data.loc[i,"ct_mean"]) > fold * data.loc[i,"ct_sd"] else True

data = data.loc[data["keep"]==True,:]

# dct
for i in sample_list:
    select = (data["samples"]==i)
    ref_mean = data.loc[select & (data["gene"]==ref),"ct_mean"].values[0]
    data.loc[select, "dct"] = data.loc[select,"ct"] - ref_mean

# ddct
## con
con_dmean = data.loc[(data["samples"]==con)].groupby("gene")["dct"].mean()

for j in gene_list:
    data.loc[(data["gene"]==j), "ddct"] = data.loc[(data["gene"]==j), "dct"] - con_dmean[j] 

# 2 ^ -ddct
data["exp"] = np.power(2, -data["ddct"])

# exp_summ
for i in sample_list:
    for j in gene_list:
        select = (data["samples"]==i) & (data["gene"] == j)
        data.loc[select,"exp_mean"] = data.loc[select,"exp"].mean()
        data.loc[select,"exp_sd"] = data.loc[select,"exp"].std()

summ = data[["samples","gene","exp_mean","exp_sd"]].drop_duplicates()
summ["upper"] = summ["exp_mean"] + summ["exp_sd"]
summ["lower"] = summ["exp_mean"] - summ["exp_sd"]

data.to_csv(outfile+".csv")

print("Plotting...")
my_color = ["#96c37d","#f3d266","#d8383a","#a9b8c6",
            "#2f7fc1","#c497b2","#b1ce46","#934b43",
            "#9394e7","#f1d77e","#d76364","#b883d3"]
# plot
p1 = (ggplot()
 +geom_bar(data=summ,mapping=aes(x="samples",y="exp_mean",fill="samples",color="samples"),alpha=0.8,stat = "identity")
 +geom_errorbar(data=summ,mapping=aes(x="samples",ymin="lower",ymax="upper",
                                        color="samples"),width=0.3)
 +geom_jitter(data=data,mapping = aes(x="samples",y="exp"),width = 0.3)
 +facet_wrap("~gene",scales = "free_y",nrow = 1)
 +theme_bw()
 +theme(panel_grid = element_blank(),
        legend_title = element_blank(),
        axis_text_x = element_text(angle = 45, hjust=1),
        strip_background = element_rect(fill = "#eaeae0"))
 +labs(x="",y="Relative expression")
 +scale_color_manual(values=my_color)
 +scale_fill_manual(values=my_color))
ggsave(self=p1,filename="Pyplot_"+outfile+".png",format="png",width=2.5*len(gene_list),
       height=5, dpi=300)
ggsave(self=p1,filename="Pyplot_"+outfile+".pdf",format="pdf",width=2.5*len(gene_list),
       height=5)
print("Done! :)")