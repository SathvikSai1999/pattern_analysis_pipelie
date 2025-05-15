import pandas as pd
from scipy.spatial import Delaunay
from scipy.sparse import csr_matrix
import csv
import numpy as np
import ast
import os

#This function is to extract the gene expression matrix of the specific pattern and the labels for DEGs analysis. 
#The output will be original pattern regions and 3-hop extended regions
def extract_motifs_id(file, gene_file, id_file,name):

    def iterate_hop(flag, iterate, mtf_simplices, tri):
        if iterate == flag:
            return
        for it_tri_sim in tri.simplices:
            #sim_name = [df['NAME'][tmp] for tmp in it_tri_sim]
            if mtf_simplices in it_tri_sim:
                for j in range(3):
                    if it_tri_sim[j] not in mtf_hop_list:
                        mtf_hop_list.append(it_tri_sim[j])
                        iterate_hop(flag, iterate + 1, it_tri_sim[j], tri)

    id_df = pd.read_csv(id_file,header=None)

    df_gene = pd.read_csv(gene_file,skiprows=1)
    gene = [row[0].split(' ') for row in df_gene.values]
    
    df = pd.read_csv(file, index_col=0)
    points = np.stack((df['X-scaled'],df['Y-scaled']),axis=-1)
    tri = Delaunay(points)
    
    ################ the index of delaunay simplices started as 0, but the df['NAME'] in the spatial file as well as the gene matrix file, the index starts from 1
    tri.simplices = [list(map(lambda x: x + 1, sim)) for sim in tri.simplices]
    
    #orginal
    allcell_label = {}
    #for i in df['NAME']:
    for i in df.index.tolist():
         allcell_label[i] = 0
    for i in range(len(id_df)):
        for j in ast.literal_eval(id_df[0][i]):
             allcell_label[j] = 1
    
    mtf_list=[i[0] for i in allcell_label.items() if i[1]==1]

    mtf_gene = [g for g in gene if int(g[1]) in mtf_list]
    #row_num = len(np.unique([row[0] for row in mtf_gene]))
    row_num = len(np.unique([row[0] for row in gene]))
    col_num = len(np.unique([row[1] for row in mtf_gene]))
    
    mtf_ids = sorted(np.unique([int(row[1]) for row in mtf_gene]))

    # Debugging output to check indices
    print("Available indices in df (first check):", df.index.tolist()[:10])  # Print first 10 row IDs
    print("Trying to filter using (first check):", mtf_ids[:10])  # Print first 10 lookup values

    # Check which IDs are missing from df.index
    missing_ids = [i for i in mtf_ids if i not in df.index]
    if missing_ids:
        print("Warning: These IDs are missing from df.index:", missing_ids)

    # Use reindex() to prevent KeyError
    filtered_df = df.reindex(mtf_ids).dropna()
    cell_types = filtered_df['top_level_cell_type']

    # Create output directory structure
    output_path = f'output/motifs/original/{name}/'
    if not os.path.exists(output_path):
         os.makedirs(output_path)
    
    with open(f'{output_path}cell_types.txt', 'w') as f:
        for item in cell_types:
            f.write(f"{item}\n")

    #reindex cell ids
    cell_reindex = {str(k):v for v,k in zip(range(1,col_num+1),mtf_ids)}
    for i in range(len(mtf_gene)):
         mtf_gene[i][1] = str(cell_reindex[mtf_gene[i][1]])
    
    motif_gene_3_cols = [' '.join(i) for i in mtf_gene]
    nonzero_num = len(mtf_gene)

    all_cell_id = sorted(allcell_label.keys())
    all_label = [[allcell_label[i]] for i in all_cell_id]

    with open(f'{output_path}mtf&non-motif_label.txt','w') as f:
         writer = csv.writer(f)
         writer.writerows(all_label)
         
    with open(f'{output_path}matrix.txt', 'w', newline='') as f:       
        writer = csv.writer(f)
        writer.writerow(['%%MatrixMarket matrix coordinate real general'])
        writer.writerow([f'{row_num} {col_num} {nonzero_num}'])
        for i in motif_gene_3_cols:
            writer.writerow([i])
    
    # 3hop
    mtf_tri = []
    
    motif_list=[i[0] for i in allcell_label.items() if i[1]==1]

    for sim in tri.simplices:
         if len(set(sim)&set(motif_list))>1:
              mtf_tri.append(sim)

    mtf_hop_list = []
    for i in motif_list:
            mtf_hop_list.append(i)

    for i,mtf in enumerate(motif_list):
            iterate_hop(3, 0, mtf, tri)

            
    for i in mtf_hop_list:
            allcell_label[i]=1

    if name.startswith('non'):
        allcell_label = {k:0 if i==1 else 1 for k,i in allcell_label.items()}
        mtf_hop_list = [k for k,i in allcell_label.items() if i==1]


    df_gene = pd.read_csv(gene_file,skiprows=1)
    gene = [row[0].split(' ') for row in df_gene.values]
    mtf_gene = [g for g in gene if int(g[1]) in mtf_hop_list]    
    row_num = len(np.unique([row[0] for row in mtf_gene]))
    col_num = len(np.unique([row[1] for row in mtf_gene]))
    
    
    mtf_ids = sorted(np.unique([int(row[1]) for row in mtf_gene]))

    # Debugging output to check indices
    print("Available indices in df (second check):", df.index.tolist()[:10])  # Print first 10 row IDs
    print("Trying to filter using (second check):", mtf_ids[:10])  # Print first 10 lookup values

    # Check which IDs are missing from df.index
    missing_ids = [i for i in mtf_ids if i not in df.index]
    if missing_ids:
        print("Warning: These IDs are missing from df.index:", missing_ids)

    # Use reindex() to prevent KeyError
    filtered_df = df.reindex(mtf_ids).dropna()
    cell_types = filtered_df['top_level_cell_type']

    # Create output directory structure for 3hop
    output_path = f'output/motifs/3hop/{name}/'
    if not os.path.exists(output_path):
         os.makedirs(output_path)
    
    with open(f'{output_path}cell_types.txt', 'w') as f:
        for item in cell_types:
            f.write(f"{item}\n")
    
    #reindex cell ids
    cell_reindex = {str(k):v for v,k in zip(range(1,col_num+1),mtf_ids)}
    for i in range(len(mtf_gene)):
         mtf_gene[i][1] = str(cell_reindex[mtf_gene[i][1]])

    motif_gene_3_cols = [' '.join(i) for i in mtf_gene]
    nonzero_num = len(mtf_gene)

    all_cell_id = sorted(allcell_label.keys())
    all_label = [[allcell_label[i]] for i in all_cell_id]
         
    with open(f'{output_path}mtf&non-motif_label.txt','w') as f:
         writer = csv.writer(f)
         writer.writerows(all_label)
         
    with open(f'{output_path}matrix.txt', 'w', newline='') as f:       
        writer = csv.writer(f)
        writer.writerow(['%%MatrixMarket matrix coordinate real general'])
        writer.writerow([f'{row_num} {col_num} {nonzero_num}'])
        for i in motif_gene_3_cols:
            writer.writerow([i])
    return


spatial_file = f'data/spatial.csv'
gene_file = f'data/matrix_raw.txt'
id_file = f'data/motif_ids.csv'

#change the motif name based on the motif you found
mtf = 'cccm'

extract_motifs_id(spatial_file, gene_file, id_file, mtf)
#change_idx(m,mtf)

