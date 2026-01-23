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
import matplotlib.pyplot as plt


# Set seeds for reproducibility
def set_seed(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
    random.seed(seed)
    torch.backends.cudnn.deterministic: bool = True
    torch.backends.cudnn.benchmark: bool = False


# Call the seed setting function
def main() -> None:
    global adata, bulk, survival_df
    
    # * set seed
    seed: int = 0
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.deterministic: bool = True
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(1)
    ite: int = 0

    set_seed(seed)

    sdh = sidish(adata, bulk, "cuda", seed=ite)