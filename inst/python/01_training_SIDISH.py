# -----------
# Initialize SIDISH with appropriate model architecture,
# learning rates, and optimizer settings.
# -------

import scanpy as sc
import pandas as pd
import numpy as np
import torch
import random
import os
import sidish


# Set seeds for reproducibility
def set_seed(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
    random.seed(seed)
    torch.backends.cudnn.deterministic: bool = True
    torch.backends.cudnn.benchmark: bool = False
    os.environ["PYTHONHASHSEED"] = str(seed)


# Call the seed setting function
def main() -> None:
    # ? get them from R scirpt
    adata: sc.AnnData = globals().get("adata")
    bulk: pd.DataFrame = globals().get("bulk")
    survival_df: pd.DataFrame = globals().get("survival_df")
    seed: int = globals().get("seed")

    # * set seed
    set_seed(seed=seed)

    ite: int = 0

    sdh = sidish(adata, bulk, "cuda", seed=ite)
