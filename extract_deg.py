import os 
import pandas as pd
import csv
import os

for region in ['3hop','original']:
    path = f'output/DEG/{region}/'
    mtf = 'cccm'
    files = f'{path}{mtf}.csv'
    df = pd.read_csv(files)
    gene_list = []
    for index, row in df.iterrows():
        if row['padj']!= 'NA':
            if float(row['padj']<0.05):
                gene_list.append([row.iloc[0]]) # Used .iloc for positional indexing


    result_path = f'output/deg_overexpressed/{region}/'
    if not os.path.exists(result_path):
         os.makedirs(result_path)
    with open(result_path+mtf+'.txt', 'w') as ff:
        writer = csv.writer(ff)
        writer.writerows(gene_list)

    
