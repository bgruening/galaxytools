import json
import os

import numpy as np
import torch
from safetensors.torch import load_file


class MinimalDataset:
    """
    Minimal dataset-like object that provides the attributes
    that Flexynesis model constructors expect.
    """

    def __init__(self, config, artifacts):
        # From config.json
        self.layers = config.get("layers", [])
        target_vars = config.get("target_variables", [])

        # Create features dict from artifacts: {layer_name: [feature1, feature2, ...]}
        self.features = artifacts["feature_lists"]

        # Create dat dict: {layer_name: tensor_placeholder}
        # Models check dat.keys() to get layer names
        self.dat = {layer: None for layer in self.layers}

        # Variable types: Variables in label_encoders are categorical,
        # variables not in label_encoders are numerical (regression)
        self.variable_types = {}

        # Annotations (ann): minimal placeholder DataFrame-like
        # Models may check np.unique(self.ann[var]) for number of classes
        self.ann = {}

        # Get class categories from label encoders in artifacts
        if "label_encoders" in artifacts:
            for var, encoder_info in artifacts["label_encoders"].items():
                if encoder_info is None:
                    continue
                if "categories" in encoder_info and len(encoder_info["categories"]) > 0:
                    # Categories are stored as [[class1, class2, ...]]
                    categories = encoder_info["categories"][0]
                    self.ann[var] = categories
                    self.variable_types[var] = "categorical"

        # For variables not in label_encoders, they are numerical (regression)
        for var in target_vars:
            if var not in self.variable_types:
                self.variable_types[var] = "numerical"
                # For numerical variables, ann should be empty or a dummy array
                # We'll use a single dummy value to indicate 1 output dimension
                self.ann[var] = np.array([0.0])


def _infer_layers_and_input_dims(config, artifacts):
    """
    Fill missing layer metadata in config using artifacts.

    Some unsupervised exports omit `layers` and `input_dims` from config.
    In that case, reconstruct them from artifacts so model modules can be
    instantiated with the same shape as the checkpoint.
    """
    feature_lists = artifacts.get("feature_lists") or {}
    if not feature_lists:
        raise ValueError(
            "artifacts.json is missing required key 'feature_lists' or it is empty"
        )

    layers = config.get("layers")
    if not layers:
        # Prefer original_modalities when present, then data_types, then feature_lists keys.
        for candidate in (
            artifacts.get("original_modalities"),
            artifacts.get("data_types"),
            list(feature_lists.keys()),
        ):
            if candidate:
                layers = list(candidate)
                break
        config["layers"] = layers
        print(f"      Inferred layers from artifacts: {layers}")

    if not layers:
        raise ValueError("Unable to infer model layers from config/artifacts")

    # Keep only layers that exist in feature_lists to avoid mismatched names.
    valid_layers = [layer for layer in layers if layer in feature_lists]
    if not valid_layers:
        raise ValueError(
            "None of the inferred/config layers exist in artifacts['feature_lists']"
        )
    if valid_layers != layers:
        print(
            "      Warning: Some layers were missing in artifacts['feature_lists']; "
            f"using {valid_layers}"
        )
        config["layers"] = valid_layers
        layers = valid_layers

    input_dims = config.get("input_dims")
    if not input_dims:
        input_dims = [len(feature_lists[layer]) for layer in layers]
        config["input_dims"] = input_dims
        print(f"      Inferred input_dims from feature_lists: {input_dims}")

    if len(input_dims) != len(layers):
        raise ValueError(
            "Length mismatch between layers and input_dims: "
            f"layers={len(layers)}, input_dims={len(input_dims)}"
        )


def reconstruct_model(safetensors_path, config_path, artifacts_path):
    """
    Reconstruct a full Flexynesis model from saved components.

    Args:
        safetensors_path: Path to .safetensors file with state_dict
        config_path: Path to config.json with model architecture
        artifacts_path: Path to artifacts.json (required)

    Returns:
        model: Fully reconstructed model instance
    """

    # 1. Load config
    print(f"[1/5] Loading config from {config_path}")
    with open(config_path, "r") as f:
        config = json.load(f)

    model_class_name = config.get("model_class")

    print(f"      Model class: {model_class_name}")
    print(f"      Input dims: {config.get('input_dims')}")
    print(f"      Layers: {config.get('layers')}")

    # 2. Load artifacts
    print(f"[2/5] Loading artifacts from {artifacts_path}")
    if not os.path.exists(artifacts_path):
        raise FileNotFoundError(f"Artifacts file not found: {artifacts_path}")
    with open(artifacts_path, "r") as f:
        artifacts = json.load(f)

    # Handle model configs (e.g. unsupervised_vae) that omit layer metadata.
    _infer_layers_and_input_dims(config, artifacts)

    # 3. Import Flexynesis model class
    print(f"[3/5] Importing Flexynesis model: {model_class_name}")

    # Only Flexynesis models are supported
    if model_class_name == "DirectPred":
        from flexynesis.models.direct_pred import DirectPred

        ModelClass = DirectPred
    elif model_class_name == "GNN":
        from flexynesis.models.gnn_early import GNN

        ModelClass = GNN
    elif model_class_name == "supervised_vae":
        from flexynesis.models.supervised_vae import supervised_vae

        ModelClass = supervised_vae
    elif model_class_name == "CrossModalPred":
        from flexynesis.models.crossmodal_pred import CrossModalPred

        ModelClass = CrossModalPred
    elif model_class_name == "MultiTripletNetwork":
        from flexynesis.models.triplet_encoder import MultiTripletNetwork

        ModelClass = MultiTripletNetwork
    else:
        raise ValueError(
            f"Unknown Flexynesis model: {model_class_name}\n"
            f"Supported models: DirectPred, GNN, supervised_vae, CrossModalPred, MultiTripletNetwork"
        )

    print(f"      Successfully imported {ModelClass.__name__}")

    # 4. Create minimal dataset object
    print("[4/5] Creating minimal dataset object")
    dataset = MinimalDataset(config, artifacts)

    # 5. Instantiate model with config
    print("[5/5] Instantiating model and loading weights")

    # Extract model config (hyperparameters)
    model_config = config.get("config", {})

    # Convert string numbers to proper types if needed
    for key in ["latent_dim", "supervisor_hidden_dim", "batch_size"]:
        if key in model_config and isinstance(model_config[key], str):
            model_config[key] = int(model_config[key])

    # Instantiate model
    model = ModelClass(
        config=model_config,
        dataset=dataset,
        target_variables=config.get("target_variables", []),
        batch_variables=None,
        surv_event_var=config.get("surv_event_var"),
        surv_time_var=config.get("surv_time_var"),
        use_loss_weighting=True,
        device_type=config.get("device_type", "cpu"),
    )

    # 6. Load state dict
    print(f"      Loading weights from SafeTensors: {safetensors_path}")
    state_dict = load_file(safetensors_path)

    model.load_state_dict(state_dict)
    model.eval()

    print("\n✓ Model reconstructed successfully!")
    print(f"  Model type: {type(model).__name__}")
    print(f"  Has .transform(): {hasattr(model, 'transform')}")
    print(f"  Has .predict(): {hasattr(model, 'predict')}")

    return model


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Reconstruct a full Flexynesis model from safetensors + config"
    )
    parser.add_argument(
        "--safetensors", required=True, help="Path to .safetensors file (state_dict)"
    )
    parser.add_argument("--config", required=True, help="Path to config.json")
    parser.add_argument("--artifacts", required=True, help="Path to artifacts.json")
    parser.add_argument(
        "--output", default="full_model.pth", help="Output path for reconstructed model"
    )

    args = parser.parse_args()

    # Reconstruct
    model = reconstruct_model(
        safetensors_path=args.safetensors,
        config_path=args.config,
        artifacts_path=args.artifacts,
    )

    # Save
    print(f"\nSaving full model to: {args.output}")
    torch.save(model, args.output)
    print("Done! Safetensors succesfully converted to Pytorch")
