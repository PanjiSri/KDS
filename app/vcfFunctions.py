import allel
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import warnings
from sklearn.preprocessing import StandardScaler


def read_vcf_for_analysis(vcf_path):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            callset = allel.read_vcf(vcf_path, fields=['samples', 'calldata/GT', 'variants/CHROM', 'variants/POS'])
            
            if callset is None:
                callset = allel.read_vcf(vcf_path, fields='*', numbers={'GT': 2})
            
            if callset is None:
                raise ValueError("Tidak dapat membaca berkas VCF. Berkas mungkin rusak atau dalam format yang tidak didukung.")
            
            if 'calldata/GT' not in callset:
                raise ValueError("Berkas VCF tidak mengandung kolom genotipe (GT) yang diperlukan untuk analisis.")
            
            if 'samples' not in callset or len(callset['samples']) == 0:
                raise ValueError("Tidak ada sampel ditemukan dalam berkas VCF.")
            
            gt_shape = callset['calldata/GT'].shape
            if gt_shape[0] == 0:
                raise ValueError("Tidak ada varian ditemukan dalam berkas VCF.")
            if gt_shape[1] == 0:
                raise ValueError("Tidak ada sampel ditemukan dalam data genotipe.")
            
            return callset
            
    except Exception as e:
        if "codec" in str(e).lower():
            raise ValueError(f"Kesalahan pengkodean berkas VCF. Pastikan berkas dalam format UTF-8: {str(e)}")
        elif "not a valid VCF" in str(e):
            raise ValueError(f"Format berkas VCF tidak valid: {str(e)}")
        else:
            raise ValueError(f"Kesalahan saat membaca berkas VCF: {str(e)}")


def apply_quality_control(callset, maf_threshold=0.05, snp_missing_threshold=0.2, ind_missing_threshold=0.2):
    try:
        gt = allel.GenotypeArray(callset['calldata/GT'])
        samples_original = callset['samples'].tolist()
        snps_original_count = gt.shape[0]
        samples_original_count = gt.shape[1]
        
        print(f"Memulai QC dengan {snps_original_count} SNP dan {samples_original_count} sampel")
        
        ac = gt.count_alleles()
        is_biallelic_01 = ac.is_biallelic_01()
        gt_bi = gt[is_biallelic_01]
        
        if gt_bi.shape[0] == 0:
            is_biallelic = ac.is_biallelic()
            if is_biallelic.sum() > 0:
                gt_bi = gt[is_biallelic]
            else:
                raise ValueError("Tidak ada SNP biallelic ditemukan dalam dataset.")
        
        print(f"Setelah filter biallelic: {gt_bi.shape[0]} SNP")
        
        ac_bi = gt_bi.count_alleles()
        
        with np.errstate(divide='ignore', invalid='ignore'): 
            aaf = ac_bi.to_frequencies()[:, 1]
            aaf = np.where(np.isnan(aaf), 0, aaf) 
        
        maf_filter = (aaf > maf_threshold) & (aaf < (1 - maf_threshold))
        gt_maf = gt_bi[maf_filter]
        
        if gt_maf.shape[0] == 0:
            max_maf_val = 0
            if len(aaf[aaf < 0.5]) > 0 : 
                max_maf_val = aaf[aaf < 0.5].max()
            raise ValueError(f"Tidak ada SNP lolos filter MAF (ambang batas: {maf_threshold}). "
                           f"MAF maksimum dalam dataset (setelah filter biallelic): {max_maf_val:.3f}")
        
        print(f"Setelah filter MAF: {gt_maf.shape[0]} SNP")
        
        missing_snp_prop = gt_maf.count_missing(axis=1) / gt_maf.shape[1]
        snp_missing_filter = missing_snp_prop < snp_missing_threshold
        gt_snp_filtered = gt_maf[snp_missing_filter]
        
        if gt_snp_filtered.shape[0] == 0:
            min_missing = 1.0 
            if len(missing_snp_prop) > 0: 
                 min_missing = missing_snp_prop.min()
            raise ValueError(f"Tidak ada SNP lolos filter data hilang (ambang batas: {snp_missing_threshold}). "
                           f"Minimum data hilang SNP (setelah filter MAF): {min_missing:.3f}")
        
        print(f"Setelah filter data hilang SNP: {gt_snp_filtered.shape[0]} SNP")
        
        missing_ind_prop = gt_snp_filtered.count_missing(axis=0) / gt_snp_filtered.shape[0]
        ind_missing_filter = missing_ind_prop < ind_missing_threshold
        gt_qc = gt_snp_filtered.compress(ind_missing_filter, axis=1)
        samples_after_qc = np.array(samples_original)[ind_missing_filter].tolist()
        
        if gt_qc.shape[1] == 0:
            min_missing = 1.0 
            if len(missing_ind_prop) > 0: 
                min_missing = missing_ind_prop.min()
            raise ValueError(f"Tidak ada sampel lolos filter data hilang (ambang batas: {ind_missing_threshold}). "
                           f"Minimum data hilang sampel (setelah filter SNP): {min_missing:.3f}")
        
        print(f"Setelah filter data hilang individu: {gt_qc.shape[1]} sampel")
        
        if gt_qc.shape[0] == 0 or gt_qc.shape[1] == 0:
            raise ValueError("Tidak ada data tersisa setelah kontrol kualitas.")
        
        snps_after_qc_count = gt_qc.shape[0]
        
        gn = gt_qc.to_n_alt(fill=-1)
        
        non_missing_gn = gn[gn != -1]
        if non_missing_gn.size == 0 or (non_missing_gn.size > 0 and non_missing_gn.std() == 0) :
            raise ValueError("Tidak ada variasi genetik tersisa setelah QC (semua genotipe sama atau hilang).")

        imputer = SimpleImputer(missing_values=-1, strategy="mean")
        gn_imputed_transposed = imputer.fit_transform(gn.T) 
        
        return gn_imputed_transposed, samples_after_qc, snps_original_count, snps_after_qc_count, samples_original_count
        
    except Exception as e:
        if isinstance(e, ValueError): 
            raise
        else: 
            raise ValueError(f"Kesalahan selama kontrol kualitas: {str(e)}")


def run_pca_analysis(genotype_matrix_imputed, n_components=3):
    try:
        n_samples, n_features = genotype_matrix_imputed.shape
        
        if n_samples < 2:
            raise ValueError(f"PCA memerlukan setidaknya 2 sampel, tetapi hanya {n_samples} tersedia.")
        
        if n_features < 1:
            raise ValueError("PCA memerlukan setidaknya 1 fitur (SNP).")
        
        effective_n_components = min(n_components, n_samples - 1, n_features)
        
        if effective_n_components < 1:
             raise ValueError(f"Tidak dapat menjalankan PCA dengan {effective_n_components} komponen. Perlu setidaknya 1 komponen yang valid berdasarkan data.")

        print(f"Menjalankan PCA dengan {effective_n_components} komponen pada {n_samples} sampel dan {n_features} fitur")
        
        scaler = StandardScaler()
        genotype_matrix_scaled = scaler.fit_transform(genotype_matrix_imputed)
        
        pca_model = PCA(n_components=effective_n_components, random_state=42)
        pcs = pca_model.fit_transform(genotype_matrix_scaled)
        explained_variance_ratio = pca_model.explained_variance_ratio_
        
        if pcs.shape[0] != n_samples:
            raise ValueError("Output PCA memiliki jumlah sampel yang salah.")
        if pcs.shape[1] != effective_n_components:
            print(f"Peringatan: PCA menghasilkan {pcs.shape[1]} komponen, diminta {effective_n_components}")
        
        if np.any(np.isnan(pcs)) or np.any(np.isinf(pcs)):
            raise ValueError("PCA menghasilkan nilai tidak valid (NaN atau Inf).")
        
        return pcs, explained_variance_ratio
        
    except Exception as e:
        if isinstance(e, ValueError):
            raise
        else:
            raise ValueError(f"Kesalahan selama analisis PCA: {str(e)}")