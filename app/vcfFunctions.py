import allel
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import warnings

def read_vcf_for_analysis(vcf_path):
    """Read VCF file using scikit-allel with better error handling."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            callset = allel.read_vcf(vcf_path, fields=['samples', 'calldata/GT', 'variants/CHROM', 'variants/POS'])
            
            if callset is None:
                callset = allel.read_vcf(vcf_path, fields='*', numbers={'GT': 2})
            
            if callset is None:
                raise ValueError("Unable to read VCF file. File may be corrupted or in an unsupported format.")
            
            if 'calldata/GT' not in callset:
                raise ValueError("VCF file does not contain genotype (GT) field required for analysis.")
            
            if 'samples' not in callset or len(callset['samples']) == 0:
                raise ValueError("No samples found in VCF file.")
            
            gt_shape = callset['calldata/GT'].shape
            if gt_shape[0] == 0:
                raise ValueError("No variants found in VCF file.")
            if gt_shape[1] == 0:
                raise ValueError("No samples found in genotype data.")
            
            return callset
            
    except Exception as e:
        if "codec" in str(e).lower():
            raise ValueError(f"VCF file encoding error. Please ensure the file is in UTF-8 format: {str(e)}")
        elif "not a valid VCF" in str(e):
            raise ValueError(f"Invalid VCF file format: {str(e)}")
        else:
            raise ValueError(f"Error reading VCF file: {str(e)}")


def apply_quality_control(callset, maf_threshold=0.05, snp_missing_threshold=0.2, ind_missing_threshold=0.2):
    """Apply quality control filters with better error handling."""
    try:
        gt = allel.GenotypeArray(callset['calldata/GT'])
        samples_original = callset['samples'].tolist()
        snps_original_count = gt.shape[0]
        samples_original_count = gt.shape[1]
        
        print(f"Starting QC with {snps_original_count} SNPs and {samples_original_count} samples")
        
        ac = gt.count_alleles()
        is_biallelic_01 = ac.is_biallelic_01()
        gt_bi = gt[is_biallelic_01]
        
        if gt_bi.shape[0] == 0:
            is_biallelic = ac.is_biallelic()
            if is_biallelic.sum() > 0:
                gt_bi = gt[is_biallelic]
            else:
                raise ValueError("No biallelic SNPs found in the dataset.")
        
        print(f"After biallelic filter: {gt_bi.shape[0]} SNPs")
        
        ac_bi = gt_bi.count_alleles()
        
        with np.errstate(divide='ignore', invalid='ignore'):
            aaf = ac_bi.to_frequencies()[:, 1]
            aaf = np.where(np.isnan(aaf), 0, aaf)
        
        maf_filter = (aaf > maf_threshold) & (aaf < (1 - maf_threshold))
        gt_maf = gt_bi[maf_filter]
        
        if gt_maf.shape[0] == 0:
            lower_threshold = maf_threshold / 2
            maf_filter = (aaf > lower_threshold) & (aaf < (1 - lower_threshold))
            gt_maf = gt_bi[maf_filter]
            if gt_maf.shape[0] == 0:
                raise ValueError(f"No SNPs passed MAF filter (threshold: {maf_threshold}). "
                               f"Maximum MAF in dataset: {aaf[aaf < 0.5].max():.3f}")
        
        print(f"After MAF filter: {gt_maf.shape[0]} SNPs")
        
        missing_snp_prop = gt_maf.count_missing(axis=1) / gt_maf.shape[1]
        snp_missing_filter = missing_snp_prop < snp_missing_threshold
        gt_snp_filtered = gt_maf[snp_missing_filter]
        
        if gt_snp_filtered.shape[0] == 0:
            min_missing = missing_snp_prop.min()
            raise ValueError(f"No SNPs passed missingness filter (threshold: {snp_missing_threshold}). "
                           f"Minimum SNP missingness: {min_missing:.3f}")
        
        print(f"After SNP missingness filter: {gt_snp_filtered.shape[0]} SNPs")
        
        missing_ind_prop = gt_snp_filtered.count_missing(axis=0) / gt_snp_filtered.shape[0]
        ind_missing_filter = missing_ind_prop < ind_missing_threshold
        gt_qc = gt_snp_filtered.compress(ind_missing_filter, axis=1)
        samples_after_qc = np.array(samples_original)[ind_missing_filter].tolist()
        
        if gt_qc.shape[1] == 0:
            min_missing = missing_ind_prop.min()
            raise ValueError(f"No samples passed missingness filter (threshold: {ind_missing_threshold}). "
                           f"Minimum sample missingness: {min_missing:.3f}")
        
        print(f"After individual missingness filter: {gt_qc.shape[1]} samples")
        
        # Final check
        if gt_qc.shape[0] == 0 or gt_qc.shape[1] == 0:
            raise ValueError("No data remaining after quality control.")
        
        snps_after_qc_count = gt_qc.shape[0]
        
        # Convert to allele counts and impute missing values
        gn = gt_qc.to_n_alt(fill=-1) 
        
        if gn[gn != -1].std() == 0:
            raise ValueError("No genetic variation remaining after QC.")
        
        imputer = SimpleImputer(missing_values=-1, strategy="mean")
        gn_imputed_transposed = imputer.fit_transform(gn.T)
        
        return gn_imputed_transposed, samples_after_qc, snps_original_count, snps_after_qc_count, samples_original_count
        
    except Exception as e:
        if isinstance(e, ValueError):
            raise
        else:
            raise ValueError(f"Error during quality control: {str(e)}")


def run_pca_analysis(genotype_matrix_imputed, n_components=3):
    """Run PCA with better validation and error handling."""
    try:
        n_samples, n_features = genotype_matrix_imputed.shape
        
        # Validate input
        if n_samples < 2:
            raise ValueError(f"PCA requires at least 2 samples, but only {n_samples} provided.")
        
        if n_features < 1:
            raise ValueError("PCA requires at least 1 feature (SNP).")
        
        # Determine actual number of components
        max_components = min(n_samples - 1, n_features, n_components)
        actual_n_components = max(1, max_components)
        
        print(f"Running PCA with {actual_n_components} components on {n_samples} samples and {n_features} features")
        
        # Standardize the data (important for PCA)
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler()
        genotype_matrix_scaled = scaler.fit_transform(genotype_matrix_imputed)
        
        # Run PCA
        pca_model = PCA(n_components=actual_n_components, random_state=42)
        pcs = pca_model.fit_transform(genotype_matrix_scaled)
        explained_variance_ratio = pca_model.explained_variance_ratio_
        
        # Validate results
        if pcs.shape[0] != n_samples:
            raise ValueError("PCA output has incorrect number of samples.")
        
        if np.any(np.isnan(pcs)) or np.any(np.isinf(pcs)):
            raise ValueError("PCA produced invalid values (NaN or Inf).")
        
        return pcs, explained_variance_ratio
        
    except Exception as e:
        if isinstance(e, ValueError):
            raise
        else:
            raise ValueError(f"Error during PCA analysis: {str(e)}")


def generate_summary_insights(variance_explained, samples_original_count, samples_after_qc_count, 
                            snps_original_count, snps_after_qc_count):
    """Generate summary insights with better formatting."""
    try:
        # Calculate filtering statistics
        samples_removed = samples_original_count - samples_after_qc_count
        snps_removed = snps_original_count - snps_after_qc_count
        samples_retention_pct = (samples_after_qc_count / samples_original_count) * 100
        snps_retention_pct = (snps_after_qc_count / snps_original_count) * 100
        
        insight_text = f"### Ringkasan Kontrol Kualitas\n\n"
        insight_text += f"**Sampel**: {samples_after_qc_count} dari {samples_original_count} "
        insight_text += f"({samples_retention_pct:.1f}% dipertahankan, {samples_removed} dihapus)\n\n"
        insight_text += f"**SNP**: {snps_after_qc_count:,} dari {snps_original_count:,} "
        insight_text += f"({snps_retention_pct:.1f}% dipertahankan, {snps_removed:,} dihapus)\n\n"
        
        insight_text += f"### Hasil PCA\n\n"
        
        if len(variance_explained) > 0:
            var_pc1 = variance_explained[0] * 100
            insight_text += f"- **PC1** menjelaskan **{var_pc1:.2f}%** dari variasi genetik\n"
        
        if len(variance_explained) > 1:
            var_pc2 = variance_explained[1] * 100
            insight_text += f"- **PC2** menjelaskan **{var_pc2:.2f}%** dari variasi genetik\n"
            cumulative_var = (variance_explained[0] + variance_explained[1]) * 100
            insight_text += f"- **PC1 & PC2** secara bersama-sama menjelaskan **{cumulative_var:.2f}%** dari total variasi\n"
        
        if len(variance_explained) > 2:
            var_pc3 = variance_explained[2] * 100
            insight_text += f"- **PC3** menjelaskan **{var_pc3:.2f}%** dari variasi genetik\n"
            cumulative_var_3 = sum(variance_explained[:3]) * 100
            insight_text += f"- **3 PC pertama** menjelaskan **{cumulative_var_3:.2f}%** dari total variasi\n"
        
        # Add interpretation
        insight_text += f"\n### Interpretasi\n\n"
        
        if len(variance_explained) > 0:
            if variance_explained[0] > 0.1:  # If PC1 explains >10%
                insight_text += "Komponen utama pertama menangkap sebagian besar variasi genetik, "
                insight_text += "menunjukkan adanya struktur populasi yang kuat atau diferensiasi dalam dataset Anda.\n"
            else:
                insight_text += "Variasi genetik tersebar di beberapa komponen, "
                insight_text += "menunjukkan adanya struktur populasi yang halus atau keragaman genetik yang tinggi.\n"
        
        return insight_text
        
    except Exception as e:
        return f"Error dalam membuat penjelasan: {str(e)}"