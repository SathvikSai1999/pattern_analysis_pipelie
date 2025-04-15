import pandas as pd
import numpy as np
from collections import Counter
from itertools import combinations
import os

def generate_rank_files(month, replicate):
    # Create result_rank directory if it doesn't exist
    if not os.path.exists('result_rank'):
        os.makedirs('result_rank')
    
    # Read cell types
    with open('data/extracted_matrix/3hop/cccm/cell_types.txt', 'r') as f:
        cell_types = [line.strip() for line in f]
    
    # Count cell type frequencies
    cell_type_counts = Counter(cell_types)
    
    # Generate node ranks
    node_ranks = []
    for cell_type, count in cell_type_counts.items():
        # Calculate rank based on frequency (higher frequency = lower rank)
        rank = len(cell_type_counts) - list(cell_type_counts.values()).index(count)
        node_ranks.append({
            'node_celltype': cell_type,
            'occurrence_num': count,
            'rank': rank
        })
    
    # Create and save node rank file
    node_df = pd.DataFrame(node_ranks)
    node_file = f'result_rank/AD_{month}m_r{replicate}_control_node_rank.csv'
    node_df.to_csv(node_file, index=False)
    print(f"Created {node_file}")
    
    # Generate tri ranks
    tri_ranks = []
    # Get unique cell types
    unique_cell_types = list(cell_type_counts.keys())
    
    # Generate all possible triplets
    for triplet in combinations(unique_cell_types, 3):
        # Count occurrences of this triplet in the data
        count = 0
        for i in range(len(cell_types)-2):
            if (cell_types[i] in triplet and 
                cell_types[i+1] in triplet and 
                cell_types[i+2] in triplet):
                count += 1
        
        if count > 0:  # Only include triplets that occur
            tri_ranks.append({
                'tri_celltype': '&'.join(sorted(triplet)),
                'occurrence_num': count,
                'rank': len(tri_ranks) + 1  # Simple rank based on order
            })
    
    # Create and save tri rank file
    tri_df = pd.DataFrame(tri_ranks)
    tri_file = f'result_rank/AD_{month}m_r{replicate}_control_tri_rank.csv'
    tri_df.to_csv(tri_file, index=False)
    print(f"Created {tri_file}")

# Generate files for all combinations
months = [8, 13]
replicates = [1, 2]

for month in months:
    for replicate in replicates:
        generate_rank_files(month, replicate) 