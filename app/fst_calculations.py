import numpy as np
import pandas as pd
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
    try:
        freq1, mask1 = calculate_allele_frequencies(df, pool1)
        freq2, mask2 = calculate_allele_frequencies(df, pool2)
        
        depth1 = df[f'pool_depth_{pool1}'].values
        depth2 = df[f'pool_depth_{pool2}'].values
        
        qualifying_snps_mask = mask1 & mask2 & (depth1 >= min_depth) & (depth2 >= min_depth)
        num_qualifying_snps = qualifying_snps_mask.sum()
        
        if num_qualifying_snps < 10:
            return np.nan
        
        p1 = freq1[qualifying_snps_mask]
        p2 = freq2[qualifying_snps_mask]
        n1 = depth1[qualifying_snps_mask]
        n2 = depth2[qualifying_snps_mask]

        p_bar = (n1 * p1 + n2 * p2) / (n1 + n2)
    
        s2_1 = p1 * (1 - p1)
        s2_2 = p2 * (1 - p2)
        s2_denominator = n1 + n2 -1
        s2_denominator_is_problem = (s2_denominator <= 0) | np.isnan(s2_denominator)
        s2 = np.full_like(s2_denominator, np.nan, dtype=float)
        s2[~s2_denominator_is_problem] = (n1[~s2_denominator_is_problem] * s2_1[~s2_denominator_is_problem] + n2[~s2_denominator_is_problem] * s2_2[~s2_denominator_is_problem]) / s2_denominator[~s2_denominator_is_problem]

        nc_denominator = n1 + n2
        nc_val = n1 + n2 - (n1**2 + n2**2) / nc_denominator
        nc_minus_1 = nc_val -1
        nc_minus_1_is_problem = (nc_minus_1 <= 0) | np.isnan(nc_minus_1)

        correction_denom = (2*(n1+n2)-1)
        correction_denom_is_problem = (correction_denom <=0) | np.isnan(correction_denom)
        correction_term = np.full_like(correction_denom, np.nan, dtype=float)
        correction_term[~correction_denom_is_problem] = (1/correction_denom[~correction_denom_is_problem]) * \
            (n1[~correction_denom_is_problem]*p1[~correction_denom_is_problem]*(1-p1[~correction_denom_is_problem]) + \
             n2[~correction_denom_is_problem]*p2[~correction_denom_is_problem]*(1-p2[~correction_denom_is_problem]))

        a = np.full_like(nc_val, np.nan, dtype=float)
        a_valid_mask = ~nc_minus_1_is_problem & ~np.isnan(s2) & ~np.isnan(correction_term)
        a[a_valid_mask] = (nc_val[a_valid_mask] / nc_minus_1[a_valid_mask]) * (s2[a_valid_mask] - correction_term[a_valid_mask])
        
        b = np.full_like(correction_term, np.nan, dtype=float)
        b_valid_mask = ~np.isnan(correction_term)
        b[b_valid_mask] = correction_term[b_valid_mask]

        c = 0.5 * (p1 - p2)**2

        numerator_fst = a + b + c
        denominator_fst = a + b + c + s2
        
        fst_per_locus = np.full_like(numerator_fst, np.nan, dtype=float)

        non_nan_components_mask = ~np.isnan(numerator_fst) & ~np.isnan(denominator_fst)
        
        calc_mask = non_nan_components_mask & (denominator_fst != 0)
        fst_per_locus[calc_mask] = numerator_fst[calc_mask] / denominator_fst[calc_mask]
        
        zero_zero_mask = non_nan_components_mask & (denominator_fst == 0) & (numerator_fst == 0)
        fst_per_locus[zero_zero_mask] = 0.0
        
        zero_nonzero_mask = non_nan_components_mask & (denominator_fst == 0) & (numerator_fst != 0)
        fst_per_locus[zero_nonzero_mask] = 1.0
        
        fst_per_locus = np.clip(fst_per_locus, 0, 1)
        
        weights = n1 + n2 
        
        valid_indices_for_avg = ~np.isnan(fst_per_locus) & ~np.isnan(weights)
        actual_fst_values = fst_per_locus[valid_indices_for_avg]
        actual_weights = weights[valid_indices_for_avg]

        if len(actual_fst_values) == 0 or np.sum(actual_weights) == 0 or np.isnan(np.sum(actual_weights)):
            return np.nan

        weighted_fst = np.sum(actual_fst_values * actual_weights) / np.sum(actual_weights)
        final_fst = np.clip(weighted_fst, 0, 1)
        return final_fst

    except Exception as e_calc:
        print(f"CRITICAL ERROR in calculate_pairwise_fst for {pool1} vs {pool2}: {type(e_calc).__name__} - {str(e_calc)}")
        import traceback
        print("Full Traceback for calculate_pairwise_fst:")
        traceback.print_exc()
        raise


def create_fst_matrix(df, min_depth=10):
    print(f"DEBUG: create_fst_matrix called. min_depth={min_depth}. df shape: {df.shape}")
    pools = extract_pool_names(df)
    n_pools = len(pools)
    print(f"DEBUG: create_fst_matrix - Number of pools: {n_pools}, Pool names: {pools}")
    
    fst_matrix = np.full((n_pools, n_pools), np.nan)
    
    for i in range(n_pools):
        fst_matrix[i, i] = 0.0
        for j in range(i + 1, n_pools):
            try:
                fst = calculate_pairwise_fst(df, pools[i], pools[j], min_depth=int(min_depth))
                fst_matrix[i, j] = fst
                fst_matrix[j, i] = fst
            except Exception as e_pair:
                print(f"ERROR: create_fst_matrix - Error during FST calculation for pair {pools[i]} vs {pools[j]}: {type(e_pair).__name__} - {str(e_pair)}")
                raise
    
    result_df = pd.DataFrame(fst_matrix, index=pools, columns=pools)
    return result_df