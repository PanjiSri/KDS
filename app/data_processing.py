import pandas as pd
import numpy as np
import io
import base64
import tempfile
import os
import warnings
from cyvcf2 import VCF 

try:
    from vcfFunctions import (
        read_vcf_for_analysis,
        apply_quality_control,
        run_pca_analysis,
        generate_summary_insights
    )
except ImportError:
    from app.vcfFunctions import (
        read_vcf_for_analysis,
        apply_quality_control,
        run_pca_analysis,
        generate_summary_insights
    )

def parse_vcf_to_json_summary(contents, filename):
    """
    Parser VCF awal menggunakan cyvcf2 untuk mendapatkan ringkasan cepat.
    Ini tidak digunakan untuk analisis utama, hanya untuk status unggah.
    """
    if contents is None:
        return None, "No VCF file uploaded."
    
    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
    except Exception as e:
        return None, f"Error decoding file: {str(e)}"
    
    # Simpan ke file temporer untuk cyvcf2
    suffix = '.vcf.gz' if filename.endswith('.gz') else '.vcf'
    is_gzipped = suffix == '.vcf.gz'

    temp_file_path = None
    try:
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix=suffix) as temp_file:
            temp_file.write(decoded)
            temp_file_path = temp_file.name
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # cyvcf2 bisa membaca .vcf atau .vcf.gz secara otomatis
            vcf_reader = VCF(temp_file_path, strict_gt=False)
            samples = list(vcf_reader.samples)
            
            # Hitung total varian (bisa lambat untuk file besar, jadi kita batasi)
            total_variants = 0
            for _ in vcf_reader:
                total_variants += 1
                if total_variants > 100000:  # Batasi iterasi untuk file besar
                    total_variants = ">100,000 (estimation)"
                    break
            vcf_reader.close()

            if not samples:
                return None, f"No samples found in VCF file '{filename}'."
            
            result_data = {
                'samples_count': len(samples),
                'total_variants_summary': str(total_variants),
                'filename': filename,
                'vcf_contents_base64': content_string
            }
            return result_data, f"Samples: {len(samples)}, Variants: {total_variants}."
        
    except Exception as e:
        print(f"VCF Summary Parsing Error for {filename}: {e}")
        error_msg = str(e)
        if "BGZF" in error_msg and not is_gzipped:
            return None, f"Error: VCF file '{filename}' appears to be gzipped but does not have a .gz extension."
        elif "not a VCF file" in error_msg or "Problem parsing" in error_msg:
            return None, f"Error: File '{filename}' does not appear to be a valid VCF file."
        elif "No such file" in error_msg:
            return None, f"Error: Could not create temporary file for processing."
        return None, f"Error parsing VCF: {error_msg}"
    finally:
        if temp_file_path and os.path.exists(temp_file_path):
            try:
                os.unlink(temp_file_path)
            except:
                pass


def trigger_analysis_pipeline(vcf_contents_base64, filename, 
                              maf_thresh=0.05, snp_miss_thresh=0.2, ind_miss_thresh=0.2, n_pca_components=3):
    """
    Orkestrasi seluruh pipeline analisis dari VCF ke PCA.
    """
    try:
        decoded_vcf = base64.b64decode(vcf_contents_base64)
    except Exception as e:
        raise RuntimeError(f"Error decoding VCF file: {str(e)}")
    
    # Simpan ke file temporer untuk scikit-allel
    suffix = '.vcf.gz' if filename.endswith('.gz') else '.vcf'
    temp_vcf_path = None
    try:
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix=suffix) as temp_file:
            temp_file.write(decoded_vcf)
            temp_vcf_path = temp_file.name
        
        # 1. Baca VCF dengan scikit-allel
        try:
            callset = read_vcf_for_analysis(temp_vcf_path)
        except Exception as e:
            raise RuntimeError(f"Error reading VCF file: {str(e)}")
        
        # 2. Quality Control
        try:
            gn_imputed_T, samples_qc, snps_orig, snps_qc, samples_orig = apply_quality_control(
                callset, 
                maf_threshold=maf_thresh, 
                snp_missing_threshold=snp_miss_thresh, 
                ind_missing_threshold=ind_miss_thresh
            )
        except Exception as e:
            raise RuntimeError(f"Error during quality control: {str(e)}")
        
        # Validate QC results
        if len(samples_qc) < 2:
            raise ValueError(f"Not enough samples for PCA after QC. Only {len(samples_qc)} samples remaining (minimum 2 required).")
        if gn_imputed_T.shape[1] < 1:
            raise ValueError(f"No SNPs remaining after QC. Try adjusting QC thresholds.")
        
        # 3. Jalankan PCA
        # Pastikan n_pca_components tidak melebihi dimensi data
        max_components = min(len(samples_qc) - 1, gn_imputed_T.shape[1])
        actual_n_components = min(n_pca_components, max_components)
        if actual_n_components < 1:
            actual_n_components = 1
        
        try:
            pcs, var_ratio = run_pca_analysis(gn_imputed_T, n_components=actual_n_components)
        except Exception as e:
            raise RuntimeError(f"Error during PCA analysis: {str(e)}")
        
        # 4. Siapkan hasil PCA untuk plot
        num_pcs_to_df = min(pcs.shape[1], n_pca_components)
        pca_columns = [f'PC{i+1}' for i in range(num_pcs_to_df)]
        df_pca_coords = pd.DataFrame(data=pcs[:, :num_pcs_to_df], columns=pca_columns)
        df_pca_coords['Sample'] = samples_qc[:len(df_pca_coords)]

        # 5. Hasilkan insight
        insights = generate_summary_insights(var_ratio, samples_orig, len(samples_qc), snps_orig, snps_qc)
        
        analysis_summary = {
            'samples_original': samples_orig,
            'samples_after_qc': len(samples_qc),
            'snps_original': snps_orig,
            'snps_after_qc': snps_qc
        }

        return {
            'pca_coords_df_json': df_pca_coords.to_json(orient='split'),
            'variance_explained': var_ratio.tolist(),
            'insights_text': insights,
            'analysis_summary': analysis_summary
        }

    except Exception as e:
        # Re-raise with more context
        if isinstance(e, (ValueError, RuntimeError)):
            raise e
        else:
            raise RuntimeError(f"Unexpected error in analysis pipeline: {str(e)}")
    finally:
        if temp_vcf_path and os.path.exists(temp_vcf_path):
            try:
                os.unlink(temp_vcf_path)
            except:
                pass


def parse_dataframe_to_json(contents, filename, file_type="CSV/TSV"):
    """Parse various file formats to JSON for Dash storage."""
    if contents is None:
        return None, f"No {file_type} file uploaded."
    
    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
    except Exception as e:
        return None, f"Error decoding file: {str(e)}"
    
    try:
        if filename.lower().endswith('.q'):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=r'\s+', header=None, engine='python')
        elif filename.lower().endswith(('.evec', '.eigen', '.pca')):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=r'\s+', header=None, engine='python')
            # Coba deteksi header atau kolom sampel
            try:
                pd.to_numeric(df.iloc[:, 0])
            except ValueError:
                df.index = df.iloc[:, 0]
                df = df.iloc[:, 1:]
        else:  # CSV atau TSV
            try:
                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            except:
                try:
                    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t')
                except:
                    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=r'\s+', engine='python')

        if df.empty:
            return None, f"{file_type} file '{filename}' is empty."

        # Validasi spesifik untuk PCA
        if file_type == "PCA":
            if df.shape[1] < 2:
                return None, f"PCA file should have at least 2 columns for PCs. Found {df.shape[1]} columns."
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) < 2:
                if df.iloc[:, 0].dtype == 'object' and df.iloc[:, 1:].select_dtypes(include=[np.number]).shape[1] >= 2:
                    pass
                else:
                    return None, f"PCA file should contain primarily numeric data for principal components."
                
        elif file_type == "ADMIXTURE":
            if df.select_dtypes(include=[np.number]).shape[1] > 0:
                numeric_data = df.select_dtypes(include=[np.number])
                if not ((numeric_data >= 0) & (numeric_data <= 1)).all().all():
                    return None, f"ADMIXTURE file should contain proportions between 0 and 1."
                row_sums = numeric_data.sum(axis=1)
                if not np.allclose(row_sums, 1.0, atol=0.01):
                    return None, f"ADMIXTURE proportions should sum to 1 for each sample."
        
        return df.to_json(date_format='iso', orient='split'), f"{file_type} file '{filename}' loaded successfully. Shape: {df.shape}."
    
    except Exception as e:
        print(f"{file_type} Parsing Error for {filename}: {e}")
        return None, f"Error parsing {file_type} file '{filename}': {str(e)}"


def parse_pca_to_json(contents, filename):
    return parse_dataframe_to_json(contents, filename, file_type="PCA")


def parse_admixture_to_json(contents, filename):
    return parse_dataframe_to_json(contents, filename, file_type="ADMIXTURE")


def parse_metadata_to_json(contents, filename):
    return parse_dataframe_to_json(contents, filename, file_type="Metadata")