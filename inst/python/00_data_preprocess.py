# -------
# 01_sidish.py
# This script receives the data from `01_sidish.R` and then executes SIDISH.
# -------

from SIDISH.SIDISH import preprocess

# * activate binding from `01_sidish.R`
global adata, bulk, survival  # ? adata$X is counts ?




def main() -> None:
    """
    Main function to preprocess the data for SIDISH.
    """
    adata, bulk_merged = preprocess(adata,
                                    bulk,
                                    survival,
                                    patient_id="Sample",
                                    celltype_name="celltype_major",
                                    processed=True)
