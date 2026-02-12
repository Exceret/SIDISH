import scanpy as sc
import pandas as pd
import numpy as np
import torch
import random
import os
from SIDISH import SIDISH as sidish
from SIDISH.SIDISH import preprocess


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
def main(
    adata: sc.AnnData = globals().get("adata"),
    bulk: pd.DataFrame = globals().get("bulk"),
    survival_df: pd.DataFrame = globals().get("survival_df"),
    seed: int = globals().get("seed"),
    save_model_dir: str = globals().get("save_model_dir")
) -> None:
    # ? get them from R scirpt

    # * set seed
    set_seed(seed=seed)

    ite: int = 0

    adata, bulk_merged = preprocess(
        adata,
        bulk,
        survival_df,
        patient_id="Sample",
        celltype_name="celltype_major",
        processed=True,
    )

    sdh: sidish = sidish(
        adata=adata,
        bulk=bulk_merged,
        device="cuda",
        seed=ite,
        use_spatial_graph=False,
        k_neighbors=None
    )

    # Reduced batch size to avoid memory issues and ensure num_workers=0

    sdh.init_Phase1(
        # Number of epochs for initial VAE training.
        epochs=225,
        # Number of iterations for retraining VAE.
        i_epochs=20,
        # Latent dimension size.
        latent_size=32,
        # List of hidden layer dimensions.
        layer_dims=[512, 128],
        batch_size=256,
        # Optimizer for VAE training.
        optimizer="Adam",
        # Learning rate.
        lr=1.0e-4,
        # Learning rate for later iterations.
        lr_3=1e-4,
        # Dropout rate.
        dropout=0,
        # Specifies dense or normal representation (default="Normal").
        type='Normal'  # or 'Dense'
    )

    sdh.init_Phase2(
        # Number of training epochs for Deep Cox model.
        epochs=500,
        # Number of neurons in the hidden layer.
        hidden=128,
        # Learning rate for Deep Cox model.
        lr=1e-4,
        # Dropout rate for training.
        dropout=0,
        # Proportion of dataset allocated to the test split.
        test_size=0.2,
        # Number of samples per batch for bulk data.
        batch_size_bulk=256
    )

    train_adata: sc.AnnData = sdh.train(
        # Number of training iterations.
        iterations=5,
        # Threshold percentile for defining High-Risk cells.
        percentile=0.95,
        # Scaling factor for updating weights.
        steepness=30,
        # Directory for saving model checkpoints.
        path=save_model_dir,
        # Number of parallel workers (default=8).
        num_workers=0,
        # If True, displays training progress (default=True).
        show=True,
        distribution_fit="fitted"  # or "default"
    )

    return train_adata


res: sc.AnnData = main(
    adata=globals().get("adata"),
    bulk=globals().get("bulk"),
    survival_df=globals().get("survival_df"),
    seed=globals().get("seed"),
    save_model_dir=globals().get("save_model_dir")
)
