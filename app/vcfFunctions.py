import allel
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

def fetchVCF(filepath):
    callset = allel.read_vcf(filepath, fields = ["samples", "calldata/GT"])
    if callset is None or "calldata/GT" not in callset:
        raise ValueError("VCF file missing 'calldata/GT' field or could not be read.")

    return callset

def imputeVCF(callset):    
    gt = allel.GenotypeArray(callset["calldata/GT"])
    gn = gt.to_n_alt()

    imputer = SimpleImputer(missing_values = -1, strategy = "mean")
    gn_imputed = imputer.fit_transform(gn)

    return gn_imputed

def convertVCFtoPCA(filepath, output_path = "pca_output.csv"):
    # Fetch Dataset
    callset = fetchVCF(filepath)
    
    # Preprocess Data
    gn_imputed = imputeVCF(callset)

    # Convert to PCA
    pca = PCA(n_components = 2)
    pcs = pca.fit_transform(gn_imputed.T)

    # Download PCA file
    df_pca = pd.DataFrame(pcs, columns=[f"PC{i + 1}" for i in range(pcs.shape[1])])
    df_pca["Sample"] = callset["samples"]

    cols = ["Sample"] + [col for col in df_pca.columns if col != "Sample"]
    df_pca = df_pca[cols]

    df_pca.to_csv(output_path, index = False)
    print(f"PCA results saved to {output_path}")

if __name__ == "__main__":
    convertVCFtoPCA("Autosomal_with_4A_markers.vcf")