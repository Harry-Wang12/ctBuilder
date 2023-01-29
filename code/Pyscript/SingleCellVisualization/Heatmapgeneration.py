from seaborn import heatmap

import numpy as np
import seaborn as sns
import matplotlib.pylab as plt

import pandas as pd


import seaborn as sns # for data visualization
import matplotlib.pyplot as plt # for data visualization

# flight = sns.load_dataset('flights') # load flights datset from GitHub seaborn repository

# # reshape flights dataeset in proper format to create seaborn heatmap
# flights_df = flight.pivot('month', 'year', 'passengers') 
# print()















dirpath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/Scoredroutes_withtreatpaper_analysis_large/"
filename = "hsa04910.xml_hsa04064.xml.txt"

scorecols = ['nashscore1','nashscore2','nashscore3','treatscore1','treatscore2','treatscore3']
pvaluecols = ['nashpvalue1','nashpvalue2','nashpvalue3','treatpvalue1','treatpvalue2','treatpvalue3']



dfdata = pd.read_csv(dirpath+filename, delimiter = "\t",usecols=pvaluecols)
# dfdata = pd.read_csv(dirpath+filename, delimiter = "\t",usecols=scorecols)


# dfdata_df = dfdata.pivot('Route ID') 










ax = sns.heatmap(dfdata, cmap=['orange','yellow'],linewidth =0.01 ,center=0.05, xticklabels=["NASH1","NASH2","NASH3","NASH Treat1","NASH Treat2","NASH Treat3"])
# ax = sns.heatmap(dfdata, cmap='PiYG_r',center=0.0, xticklabels=["NASH1","NASH2","NASH3","NASH Treat1","NASH Treat2","NASH Treat3"])
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
ax.set_ylabel("Route ID")
ax.set_title("Crosstalk routes score pvalue heatmap")
plt.show()