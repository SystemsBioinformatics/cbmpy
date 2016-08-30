import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re

fn = 'Result_meth0-(bulgTherSrcOnlyBulg)-(consets)-(minset)-(14.0)-(128).csv'

data = pd.read_csv(fn)

dataMs = data[data['Minset'] == 1]
dataMs.drop(['Minset', 'OptFlux', 'MinSumFlux', 'Constraint'], 1, inplace=True)
dataMs.set_index('R', inplace=True)
dataMs[dataMs < 0] = -1
dataMs[dataMs > 0] = 1

dataMs.sort_index(inplace=True)
fontsize_pt = plt.rcParams['ytick.labelsize']
dpi = 72.27

# comput the matrix height in points and inches
matrix_height_pt = fontsize_pt * dataMs.shape[0]
matrix_height_in = matrix_height_pt / dpi

# compute the required figure height
top_margin = 0.04  # in percentage of the figure height
bottom_margin = 0.24 # in percentage of the figure height
figure_height = matrix_height_in / (1 - top_margin - bottom_margin)


# build the figure instance with the desired height
fig, ax = plt.subplots(
        figsize=(18, figure_height),
        gridspec_kw=dict(top=1-top_margin, bottom=bottom_margin))

# let seaborn do it's thing
ax = sns.heatmap(dataMs, ax=ax)  #, cmap='jet')
for item in ax.get_yticklabels():
    item.set_rotation(0)
for item in ax.get_xticklabels():
    item.set_rotation(90)

plt.savefig(re.sub('.csv', '.png', fn))
plt.show()
