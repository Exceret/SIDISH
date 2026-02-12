from typing import Any
from typing import Callable
import inspect
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


def filter_args_4_func(input_dict: dict[str, Any], func: Callable) -> dict[str, Any]:
    """
    从 input_dict 中提取出所有键名与 func 的参数名匹配的键值对。

    Parameters
    ----------
    input_dict : dict
        输入字典，包含可能的参数键值对。
    func : callable
        目标函数。
    Returns
    -------
    dict
        仅包含 func 所需参数的子集字典，可用于 **kwargs 调用。
    """
    sig = inspect.signature(func)
    func_params = set(sig.parameters.keys())
    return {k: v for k, v in input_dict.items() if k in func_params}


# Call the seed setting function
def main(
    # ? get them from R scirpt
    adata: sc.AnnData,
    bulk: pd.DataFrame,
    survival_df: pd.DataFrame,
    seed: int = 123,
) -> None:
    other_args: dict[str, Any] = globals()

    # * set seed
    set_seed(seed=seed)

    ite: int = 0

    adata, bulk_merged = preprocess(
        adata=adata,
        bulk=bulk,
        survival_df=survival_df,
        patient_id=other_args.get("patient_id", "Sample"),
        celltype_name=other_args.get("celltype_name", "celltype_major"),
        processed=other_args.get("processed", True),
    )

    sdh: sidish = sidish(
        adata=adata,
        bulk=bulk_merged,
        device=other_args.get("device", "cuda"),
        seed=ite,
        use_spatial_graph=other_args.get("use_spatial_graph", False),
        k_neighbors=other_args.get("k_neighbors", None),
    )

    # Reduced batch size to avoid memory issues and ensure num_workers=0

    sdh.init_Phase1(
        # Number of epochs for initial VAE training.
        epochs=other_args.get("phase1_epochs", 225),
        # Number of iterations for retraining VAE.
        i_epochs=other_args.get("phase1_i_epochs", 20),
        # Latent dimension size.
        latent_size=other_args.get("phase1_latent_size", 32),
        # List of hidden layer dimensions.
        layer_dims=other_args.get("phase1_layer_dims", [512, 128]),
        batch_size=other_args.get("phase1_batch_size", 256),
        # Optimizer for VAE training.
        optimizer=other_args.get("phase1_optimizer", "Adam"),
        # Learning rate.
        lr=other_args.get("phase1_lr", 1.0e-4),
        # Learning rate for later iterations.
        lr_3=other_args.get("phase1_lr_3", 1e-4),
        # Dropout rate.
        dropout=other_args.get("phase1_dropout", 0),
        # Specifies dense or normal representation (default="Normal").
        type=other_args.get("phase1_type", "Normal"),  # or 'Dense'
    )

    sdh.init_Phase2(
        # Number of training epochs for Deep Cox model.
        epochs=other_args.get("phase2_epochs", 500),
        # Number of neurons in the hidden layer.
        hidden=other_args.get("phase2_hidden", 128),
        # Learning rate for Deep Cox model.
        lr=other_args.get("phase2_lr", 1e-4),
        # Dropout rate for training.
        dropout=other_args.get("phase2_dropout", 0),
        # Proportion of dataset allocated to the test split.
        test_size=other_args.get("phase2_test_size", 0.2),
        # Number of samples per batch for bulk data.
        batch_size_bulk=other_args.get("phase2_batch_size_bulk", 256),
    )

    train_adata: sc.AnnData = sdh.train(
        # Number of training iterations.
        iterations=other_args.get("train_iterations", 5),
        # Threshold percentile for defining High-Risk cells.
        percentile=other_args.get("train_percentile", 0.95),
        # Scaling factor for updating weights.
        steepness=other_args.get("train_steepness", 30),
        # Directory for saving model checkpoints.
        path=other_args.get("train_path", "./"),
        # Number of parallel workers (default=8).
        num_workers=other_args.get("train_num_workers", 0),
        # If True, displays training progress (default=True).
        show=other_args.get("verbose", True),
        distribution_fit=other_args.get(
            "train_distribution_fit", "fitted"
        ),  # or "default"
    )

    return train_adata


# ? pass to R script
res: sc.AnnData = main(
    adata=globals().get("adata"),
    bulk=globals().get("bulk"),
    survival_df=globals().get("survival_df"),
    seed=globals().get("seed"),
    save_model_dir=globals().get("save_model_dir"),
)

obs: pd.DataFrame = res.obs
