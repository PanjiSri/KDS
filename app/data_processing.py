import pandas as pd
import numpy as np
import io
import base64
import tempfile
import os
import warnings
from cyvcf2 import VCF

def parse_vcf_to_json(contents, filename):
    if contents is None:
        return None, "No VCF file uploaded."
    
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    
    with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.vcf') as temp_file:
        temp_file.write(decoded)
        temp_file_path = temp_file.name
    
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*Contig.*is not defined in the header.*")
            
            vcf_reader = VCF(temp_file_path, strict_gt=False)
            samples = list(vcf_reader.samples)
            variants_data = []
            total_variants = 0
            
            genotype_preview = []
            
            for i, variant in enumerate(vcf_reader):
                total_variants += 1
                if i < 5:  
                    variants_data.append({
                        'CHROM': variant.CHROM, 
                        'POS': variant.POS, 
                        'REF': variant.REF, 
                        'ALT': ','.join(variant.ALT) if variant.ALT else ''
                    })
                
                if i < 100:
                    gts = variant.gt_types  
                    genotype_preview.append(gts.tolist())
            
            vcf_reader.close()
            
            if not samples:
                return None, f"No samples found in VCF file '{filename}'."
            
            result_data = {
                'samples': samples,
                'total_variants': total_variants,
                'variants_preview': variants_data,
                'genotype_preview': genotype_preview
            }
            
            return result_data, f"VCF file '{filename}' processed successfully. Samples: {len(samples)}, Total variants: {total_variants}."
        
    except Exception as e:
        print(f"VCF Parsing Error for {filename}: {e}")
        if "BGZF" in str(e):
            return None, f"Error: VCF file appears to be compressed. Please upload an uncompressed VCF file."
        elif "truncated" in str(e):
            return None, f"Error: VCF file appears to be truncated or corrupted."
        else:
            return None, f"Error parsing VCF file '{filename}': {str(e)}"
    
    finally:
        try:
            os.unlink(temp_file_path)
        except:
            pass

def parse_dataframe_to_json(contents, filename, file_type="CSV/TSV"):
    if contents is None:
        return None, f"No {file_type} file uploaded."
    
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    
    try:
        if filename.lower().endswith('.q'):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=r'\s+', header=None, engine='python')
        elif filename.lower().endswith(('.evec', '.eigen')):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=r'\s+', header=None, engine='python')

            try:
                pd.to_numeric(df.iloc[:, 0])
            except:
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
            return None, f"{file_type} file '{filename}' is empty."

        if file_type == "PCA":
            if df.shape[1] < 2:
                return None, f"PCA file should have at least 2 columns (PC1, PC2). Found {df.shape[1]} columns."
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) < 2:
                return None, f"PCA file should contain numeric data for principal components."
                
        elif file_type == "ADMIXTURE":
            if df.select_dtypes(include=[np.number]).shape[1] > 0:
                numeric_data = df.select_dtypes(include=[np.number])
                if not ((numeric_data >= 0) & (numeric_data <= 1)).all().all():
                    return None, f"ADMIXTURE file should contain proportions between 0 and 1."
                row_sums = numeric_data.sum(axis=1)
                if not np.allclose(row_sums, 1.0, atol=0.01):
                    return None, f"ADMIXTURE proportions should sum to 1 for each sample. Found sums ranging from {row_sums.min():.3f} to {row_sums.max():.3f}."
        
        return df.to_json(date_format='iso', orient='split'), f"{file_type} file '{filename}' loaded successfully. Shape: {df.shape}."
    
    except Exception as e:
        print(f"{file_type} Parsing Error for {filename}: {e}")
        return None, f"Error parsing {file_type} file '{filename}': {str(e)}. Please check the file format."

def parse_pca_to_json(contents, filename):
    return parse_dataframe_to_json(contents, filename, file_type="PCA")

def parse_admixture_to_json(contents, filename):
    return parse_dataframe_to_json(contents, filename, file_type="ADMIXTURE")

def parse_metadata_to_json(contents, filename):
    return parse_dataframe_to_json(contents, filename, file_type="Metadata")