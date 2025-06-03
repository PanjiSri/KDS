import pandas as pd
import numpy as np
import io
import base64
import tempfile
import os
import warnings
from cyvcf2 import VCF
import traceback

try:
    from vcfFunctions import (
        read_vcf_for_analysis,
        apply_quality_control,
        run_pca_analysis
    )
    from fst_calculations import (
        create_fst_matrix
    )
except ImportError:
    from app.vcfFunctions import (
        read_vcf_for_analysis,
        apply_quality_control,
        run_pca_analysis
    )
    from app.fst_calculations import (
        create_fst_matrix
    )

def parse_vcf_to_json_summary(contents, filename):
    if contents is None:
        return None, "Tidak ada berkas VCF yang diunggah."
    
    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
    except Exception as e:
        return None, f"Kesalahan saat mendekode berkas: {str(e)}"
    
    suffix = '.vcf.gz' if filename.endswith('.gz') else '.vcf'
    is_gzipped = suffix == '.vcf.gz'

    temp_file_path = None
    try:
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix=suffix) as temp_file:
            temp_file.write(decoded)
            temp_file_path = temp_file.name
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            vcf_reader = VCF(temp_file_path, strict_gt=False)
            samples = list(vcf_reader.samples)
            
            total_variants = 0
            for _ in vcf_reader:
                total_variants += 1
                if total_variants > 100000:
                    total_variants = ">100.000 (estimasi)"
                    break
            vcf_reader.close()

            if not samples:
                return None, f"Tidak ada sampel ditemukan dalam berkas VCF '{filename}'."
            
            result_data = {
                'samples_count': len(samples),
                'total_variants_summary': str(total_variants),
                'filename': filename,
                'vcf_contents_base64': content_string
            }
            return result_data, f"Sampel: {len(samples)}, Varian: {total_variants}."
        
    except Exception as e:
        error_msg = str(e)
        if "BGZF" in error_msg and not is_gzipped:
            return None, f"Kesalahan: Berkas VCF '{filename}' tampak terkompresi gzip tetapi tidak memiliki ekstensi .gz."
        elif "not a VCF file" in error_msg or "Problem parsing" in error_msg:
            return None, f"Kesalahan: Berkas '{filename}' bukan merupakan berkas VCF yang valid."
        elif "No such file" in error_msg:
            return None, f"Kesalahan: Tidak dapat membuat berkas temporer untuk pemrosesan."
        return None, f"Kesalahan saat memparsing VCF: {error_msg}"
    finally:
        if temp_file_path and os.path.exists(temp_file_path):
            try:
                os.unlink(temp_file_path)
            except:
                pass


def trigger_analysis_pipeline(vcf_contents_base64, filename,
                              maf_thresh=0.05, snp_miss_thresh=0.2, ind_miss_thresh=0.2, n_pca_components=3):
    try:
        decoded_vcf = base64.b64decode(vcf_contents_base64)
    except Exception as e:
        raise RuntimeError(f"Kesalahan saat mendekode berkas VCF: {str(e)}")
    
    suffix = '.vcf.gz' if filename.endswith('.gz') else '.vcf'
    temp_vcf_path = None
    try:
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix=suffix) as temp_file:
            temp_file.write(decoded_vcf)
            temp_vcf_path = temp_file.name
        
        try:
            callset = read_vcf_for_analysis(temp_vcf_path)
        except Exception as e:
            raise RuntimeError(f"Kesalahan saat membaca berkas VCF: {str(e)}")
        
        try:
            gn_imputed_T, samples_qc, snps_orig, snps_qc, samples_orig = apply_quality_control(
                callset,
                maf_threshold=maf_thresh,
                snp_missing_threshold=snp_miss_thresh,
                ind_missing_threshold=ind_miss_thresh
            )
        except Exception as e:
            raise RuntimeError(f"Kesalahan selama kontrol kualitas: {str(e)}")
        
        if len(samples_qc) < 2: 
            raise ValueError(f"Jumlah sampel tidak cukup untuk PCA setelah QC. Hanya {len(samples_qc)} sampel tersisa (minimal 2 dibutuhkan).")
        if gn_imputed_T.shape[1] < 1: 
            raise ValueError(f"Tidak ada SNP tersisa setelah QC. Coba sesuaikan ambang batas QC.")
        
        max_possible_components = min(len(samples_qc) - 1, gn_imputed_T.shape[1])
        
        actual_n_components = min(n_pca_components, max_possible_components)
        
        if actual_n_components < 1:
            if max_possible_components >=1:
                actual_n_components = 1
            else: 
                 raise ValueError(f"Data tidak cukup untuk menghitung komponen PCA (sampel setelah QC={len(samples_qc)}, fitur setelah QC={gn_imputed_T.shape[1]}). Perlu minimal 2 sampel dan 1 fitur.")

        try:
            pcs, var_ratio = run_pca_analysis(gn_imputed_T, n_components=actual_n_components)
        except Exception as e:
            raise RuntimeError(f"Kesalahan selama analisis PCA: {str(e)}")
        
        num_pcs_generated = pcs.shape[1] 
        pca_columns = [f'PC{i+1}' for i in range(num_pcs_generated)]
        df_pca_coords = pd.DataFrame(data=pcs[:, :num_pcs_generated], columns=pca_columns)
        df_pca_coords['Sample'] = samples_qc[:len(df_pca_coords)]
        
        analysis_summary = {
            'samples_original': samples_orig,
            'samples_after_qc': len(samples_qc),
            'snps_original': snps_orig,
            'snps_after_qc': snps_qc
        }

        return {
            'pca_coords_df_json': df_pca_coords.to_json(orient='split'),
            'variance_explained': var_ratio.tolist(),
            'analysis_summary': analysis_summary
        }

    except Exception as e:
        if isinstance(e, (ValueError, RuntimeError)):
            raise e
        else:
            raise RuntimeError(f"Kesalahan tak terduga dalam pipeline analisis: {str(e)}")
    finally:
        if temp_vcf_path and os.path.exists(temp_vcf_path):
            try:
                os.unlink(temp_vcf_path)
            except:
                pass


def parse_dataframe_to_json(contents, filename, file_type="CSV/TSV"):
    if contents is None:
        return None, f"Tidak ada berkas {file_type} yang diunggah."
    
    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
    except Exception as e:
        return None, f"Kesalahan saat mendekode berkas: {str(e)}"
    
    try:
        if filename.lower().endswith('.q'):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=r'\s+', header=None, engine='python')
        elif filename.lower().endswith(('.evec', '.eigen', '.pca')):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=r'\s+', header=None, engine='python')
            try:
                pd.to_numeric(df.iloc[:, 0])
            except ValueError:
                df.index = df.iloc[:, 0]
                df = df.iloc[:, 1:]
        else:
            try:
                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            except:
                try:
                    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t')
                except:
                    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=r'\s+', engine='python')

        if df.empty:
            return None, f"Berkas {file_type} '{filename}' kosong."

        if file_type == "PCA":
            if df.shape[1] < 2:
                return None, f"Berkas PCA harus memiliki setidaknya 2 kolom (Sampel, PC1). Ditemukan {df.shape[1]} kolom."
            numeric_cols = df.iloc[:,1:].select_dtypes(include=[np.number]).columns
            if len(numeric_cols) < 1: 
                 return None, f"Berkas PCA harus berisi data numerik untuk komponen utama."
                
        elif file_type == "ADMIXTURE":
            if df.select_dtypes(include=[np.number]).shape[1] > 0:
                numeric_data = df.select_dtypes(include=[np.number])
                if not ((numeric_data >= 0) & (numeric_data <= 1)).all().all():
                    return None, f"Berkas ADMIXTURE harus berisi proporsi antara 0 dan 1."
                row_sums = numeric_data.sum(axis=1)
                if not np.allclose(row_sums, 1.0, atol=0.01): 
                    return None, f"Proporsi ADMIXTURE harus berjumlah mendekati 1 untuk setiap sampel."
        
        return df.to_json(date_format='iso', orient='split'), f"Berkas {file_type} '{filename}' berhasil dimuat. Bentuk: {df.shape}."
    
    except Exception as e:
        return None, f"Kesalahan saat memparsing berkas {file_type} '{filename}': {str(e)}"


def parse_pca_to_json(contents, filename):
    return parse_dataframe_to_json(contents, filename, file_type="PCA")


def parse_admixture_to_json(contents, filename):
    return parse_dataframe_to_json(contents, filename, file_type="ADMIXTURE")


def parse_metadata_to_json(contents, filename):
    return parse_dataframe_to_json(contents, filename, file_type="Metadata")


def parse_pooled_data(contents, filename):
    if contents is None:
        return None, "Tidak ada berkas pooled data yang diunggah."
    
    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
    except Exception as e:
        return None, f"Kesalahan saat mendekode berkas: {str(e)}"
    
    try:
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=r'\s+', engine='python')

        if df.empty:
            return None, "Berkas pooled data kosong."
        
        all_cols = df.columns.tolist()
        ref_cols = [col for col in all_cols if col.startswith('reference_count_')]
        depth_cols = [col for col in all_cols if col.startswith('pool_depth_')]
        
        if not ref_cols or not depth_cols:
            return None, "Berkas tidak memiliki kolom reference_count_ atau pool_depth_ yang diperlukan."

        critical_cols = ref_cols + depth_cols
        for col_name in critical_cols:
            if col_name in df.columns:
                original_dtype = df[col_name].dtype
                df[col_name] = pd.to_numeric(df[col_name], errors='coerce')
                if df[col_name].isnull().any():
                    print(f"WARNING: parse_pooled_data - Column '{col_name}' (original dtype: {original_dtype}) contains NaNs after pd.to_numeric(errors='coerce'). Check input file for non-numeric values in this column.")
            else:
                return None, f"Kolom kritis '{col_name}' tidak ditemukan."

        for col_name in critical_cols:
            if not pd.api.types.is_numeric_dtype(df[col_name]):
                 return None, f"Kolom kritis {col_name} gagal dikonversi ke tipe numerik."

        for col_name in depth_cols:
            if (df[col_name] < 0).any():
                return None, f"Kolom kedalaman '{col_name}' berisi nilai negatif yang tidak valid."

        pool_suffixes = set()
        for c in ref_cols: pool_suffixes.add(c.replace('reference_count_',''))
        for c in depth_cols: pool_suffixes.add(c.replace('pool_depth_',''))

        for suffix in sorted(list(pool_suffixes)):
            ref_col = f'reference_count_{suffix}'
            depth_col = f'pool_depth_{suffix}'
            if ref_col in df.columns and depth_col in df.columns:
                invalid_rows_mask = (df[ref_col] > df[depth_col]) & df[depth_col].notna() & df[ref_col].notna()
                if invalid_rows_mask.any():
                    return None, f"Untuk pool {suffix}, ditemukan reference_count > pool_depth pada beberapa baris. Ini tidak valid."

        return df.to_json(orient='split'), f"Berkas pooled data '{filename}' berhasil dimuat. SNPs: {len(df)}, Pools: {len(ref_cols)}."
    
    except Exception as e:
        print(f"ERROR: Unhandled exception in parse_pooled_data for {filename}: {type(e).__name__} - {str(e)}")
        print("Full Traceback for parse_pooled_data:")
        traceback.print_exc()
        return None, f"Kesalahan saat memparsing berkas pooled data: {str(e)}"


def analyze_fst_from_pooled_data(pooled_df_json, min_depth=10):
    try:
        df = pd.read_json(io.StringIO(pooled_df_json), orient='split')
        
        fst_matrix = create_fst_matrix(df, min_depth=int(min_depth))
        
        if fst_matrix is None:
             raise RuntimeError("create_fst_matrix returned None, which is unexpected.")

        return {
            'fst_matrix': fst_matrix.to_json(orient='split')
        }
    
    except Exception as e:
        print(f"ERROR: Exception in analyze_fst_from_pooled_data: {type(e).__name__} - {str(e)}")
        print("Full Traceback for analyze_fst_from_pooled_data:")
        traceback.print_exc()
        raise