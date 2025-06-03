import numpy as np
import pandas as pd
# from scipy.spatial.distance import squareform
# from scipy.cluster.hierarchy import linkage, dendrogram 
import warnings


def extract_pool_names(df):
    ref_cols = [col for col in df.columns if col.startswith('reference_count_')]
    pool_names = [col.replace('reference_count_', '') for col in ref_cols]
    return sorted(pool_names)


def calculate_allele_frequencies(df, pool_name):
    ref_col = f'reference_count_{pool_name}'
    depth_col = f'pool_depth_{pool_name}'
    
    if ref_col not in df.columns or depth_col not in df.columns:
        raise ValueError(f"Missing columns for pool {pool_name}")
    
    ref_counts = df[ref_col].values
    depths = df[depth_col].values
    
    valid_mask = depths > 0
    freqs = np.zeros(len(df))
    freqs[valid_mask] = ref_counts[valid_mask] / depths[valid_mask]
    
    return freqs, valid_mask


def calculate_pairwise_fst(df, pool1, pool2, min_depth=10):
    freq1, mask1 = calculate_allele_frequencies(df, pool1)
    freq2, mask2 = calculate_allele_frequencies(df, pool2)
    
    depth1 = df[f'pool_depth_{pool1}'].values
    depth2 = df[f'pool_depth_{pool2}'].values
    
    valid_mask = mask1 & mask2 & (depth1 >= min_depth) & (depth2 >= min_depth)
    
    if valid_mask.sum() < 10:
        return np.nan
    
    p1 = freq1[valid_mask]
    p2 = freq2[valid_mask]
    n1 = depth1[valid_mask]
    n2 = depth2[valid_mask]
    
    p_bar = (n1 * p1 + n2 * p2) / (n1 + n2)
    
    s2_1 = p1 * (1 - p1)
    s2_2 = p2 * (1 - p2)
    s2 = (n1 * s2_1 + n2 * s2_2) / (n1 + n2 - 1)
    
    nc = n1 + n2 - (n1**2 + n2**2) / (n1 + n2)
    
    a = nc / (nc - 1) * (s2 - 1/(2*(n1+n2)-1) * ((n1*p1*(1-p1) + n2*p2*(1-p2))))
    b = 1 / (2*(n1+n2)-1) * (n1*p1*(1-p1) + n2*p2*(1-p2))
    c = 0.5 * (p1 - p2)**2
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fst_per_snp = (a + b + c) / (a + b + c + s2)
        fst_per_snp = np.clip(fst_per_snp, 0, 1)
    
    weights = n1 + n2
    weighted_fst = np.sum(fst_per_snp * weights) / np.sum(weights)
    
    return np.clip(weighted_fst, 0, 1)


def create_fst_matrix(df, min_depth=10):
    pools = extract_pool_names(df)
    n_pools = len(pools)
    
    fst_matrix = np.zeros((n_pools, n_pools))
    
    for i in range(n_pools):
        for j in range(i+1, n_pools):
            fst = calculate_pairwise_fst(df, pools[i], pools[j], min_depth)
            if not np.isnan(fst):
                fst_matrix[i, j] = fst
                fst_matrix[j, i] = fst
    
    return pd.DataFrame(fst_matrix, index=pools, columns=pools)


# def calculate_dendrogram_data(fst_matrix_df):
#     pools = fst_matrix_df.index.tolist()
#     matrix_values = fst_matrix_df.values
    
#     condensed_matrix = squareform(matrix_values)
    
#     linkage_matrix = linkage(condensed_matrix, method='average')
    
#     return linkage_matrix, pools