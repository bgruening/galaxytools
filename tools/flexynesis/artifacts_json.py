"""
JSON serialization utilities for flexynesis artifacts.
Provides pickle-free alternative for Galaxy platform compatibility.

Supports conversion of sklearn objects:
- StandardScaler
- LabelEncoder
- OrdinalEncoder
"""

import json
import numpy as np
from sklearn.preprocessing import StandardScaler, LabelEncoder, OrdinalEncoder


def artifacts_to_json(artifacts, output_path):
    """
    Convert artifacts dictionary (with sklearn objects) to JSON.

    Args:
        artifacts (dict): Artifacts dictionary from training with sklearn objects
        output_path (str): Path to save JSON file

    Example:
        >>> artifacts = joblib.load('model.artifacts.joblib')
        >>> artifacts_to_json(artifacts, 'model.artifacts.json')
    """
    json_artifacts = {
        'schema_version': artifacts.get('schema_version', 1),
        'data_types': artifacts.get('data_types', []),
        'original_modalities': artifacts.get('original_modalities', []),
        'target_variables': artifacts.get('target_variables', []),
        'covariate_vars': artifacts.get('covariate_vars', []),
        'join_key': artifacts.get('join_key', None),
        'string_organism': artifacts.get('string_organism', None),
        'string_node_name': artifacts.get('string_node_name', None),
        'feature_lists': {},
        'transforms': {},
        'label_encoders': {}
    }

    # Convert feature lists (already JSON-compatible)
    feature_lists = artifacts.get('feature_lists', {})
    for modality, features in feature_lists.items():
        json_artifacts['feature_lists'][modality] = list(features)

    # Convert StandardScaler objects to JSON-compatible dicts
    transforms = artifacts.get('transforms', {})
    for modality, scaler in transforms.items():
        if scaler is None:
            json_artifacts['transforms'][modality] = None
            continue

        # Extract parameters from StandardScaler
        scaler_dict = {
            'type': 'StandardScaler',
            'with_mean': scaler.with_mean,
            'with_std': scaler.with_std,
        }

        # Convert numpy arrays to lists
        if hasattr(scaler, 'mean_'):
            scaler_dict['mean'] = scaler.mean_.tolist()
        if hasattr(scaler, 'scale_'):
            scaler_dict['scale'] = scaler.scale_.tolist()
        if hasattr(scaler, 'var_'):
            scaler_dict['var'] = scaler.var_.tolist()
        if hasattr(scaler, 'n_features_in_'):
            scaler_dict['n_features_in'] = int(scaler.n_features_in_)
        if hasattr(scaler, 'feature_names_in_'):
            scaler_dict['feature_names_in'] = scaler.feature_names_in_.tolist()
        if hasattr(scaler, 'n_samples_seen_'):
            scaler_dict['n_samples_seen'] = int(scaler.n_samples_seen_)

        json_artifacts['transforms'][modality] = scaler_dict

    # Convert LabelEncoder/OrdinalEncoder objects to JSON-compatible dicts
    label_encoders = artifacts.get('label_encoders', {})
    for variable, encoder in label_encoders.items():
        if encoder is None:
            json_artifacts['label_encoders'][variable] = None
            continue

        # Handle LabelEncoder
        if isinstance(encoder, LabelEncoder):
            encoder_dict = {
                'type': 'LabelEncoder',
                'classes': encoder.classes_.tolist()
            }
        # Handle OrdinalEncoder
        elif isinstance(encoder, OrdinalEncoder):
            encoder_dict = {
                'type': 'OrdinalEncoder',
                'categories': [cat.tolist() for cat in encoder.categories_],
            }
            if hasattr(encoder, 'handle_unknown'):
                encoder_dict['handle_unknown'] = encoder.handle_unknown
            if hasattr(encoder, 'unknown_value'):
                encoder_dict['unknown_value'] = encoder.unknown_value
            if hasattr(encoder, 'encoded_missing_value'):
                encoder_dict['encoded_missing_value'] = encoder.encoded_missing_value
            if hasattr(encoder, 'n_features_in_'):
                encoder_dict['n_features_in'] = int(encoder.n_features_in_)
            if hasattr(encoder, 'feature_names_in_'):
                encoder_dict['feature_names_in'] = encoder.feature_names_in_.tolist()
            if hasattr(encoder, '_missing_indices'):
                encoder_dict['_missing_indices'] = encoder._missing_indices
            if hasattr(encoder, '_infrequent_enabled'):
                encoder_dict['_infrequent_enabled'] = encoder._infrequent_enabled
        else:
            raise ValueError(f"Unknown encoder type: {type(encoder).__name__}")
        
        json_artifacts['label_encoders'][variable] = encoder_dict

    # Save to JSON file
    with open(output_path, 'w') as f:
        json.dump(json_artifacts, f, indent=2)

    print(f"[INFO] Saved artifacts to JSON: {output_path}")
    return json_artifacts


def json_to_artifacts(json_path):
    """
    Convert JSON artifacts back to dictionary with reconstructed sklearn objects.

    Args:
        json_path (str): Path to JSON artifacts file

    Returns:
        dict: Artifacts dictionary with reconstructed sklearn objects

    Example:
        >>> artifacts = json_to_artifacts('model.artifacts.json')
        >>> scaler = artifacts['transforms']['gex']
        >>> scaled_data = scaler.transform(test_data)
    """
    # Load JSON
    with open(json_path, 'r') as f:
        json_artifacts = json.load(f)

    artifacts = {
        'schema_version': json_artifacts.get('schema_version', 1),
        'data_types': json_artifacts.get('data_types', []),
        'original_modalities': json_artifacts.get('original_modalities', []),
        'target_variables': json_artifacts.get('target_variables', []),
        'covariate_vars': json_artifacts.get('covariate_vars', []),
        'join_key': json_artifacts.get('join_key', None),
        'string_organism': json_artifacts.get('string_organism', None),
        'string_node_name': json_artifacts.get('string_node_name', None),
        'feature_lists': {},
        'transforms': {},
        'label_encoders': {}
    }

    # Feature lists are already in correct format
    artifacts['feature_lists'] = json_artifacts.get('feature_lists', {})

    # Reconstruct StandardScaler objects
    transforms_json = json_artifacts.get('transforms', {})
    for modality, scaler_dict in transforms_json.items():
        if scaler_dict is None:
            artifacts['transforms'][modality] = None
            continue

        if scaler_dict.get('type') != 'StandardScaler':
            raise ValueError(f"Unknown scaler type: {scaler_dict.get('type')}")

        # Create StandardScaler and set its attributes
        scaler = StandardScaler(
            with_mean=scaler_dict.get('with_mean', True),
            with_std=scaler_dict.get('with_std', True)
        )

        # Restore fitted attributes
        if 'mean' in scaler_dict:
            scaler.mean_ = np.array(scaler_dict['mean'])
        if 'scale' in scaler_dict:
            scaler.scale_ = np.array(scaler_dict['scale'])
        if 'var' in scaler_dict:
            scaler.var_ = np.array(scaler_dict['var'])
        if 'n_features_in' in scaler_dict:
            scaler.n_features_in_ = scaler_dict['n_features_in']
        if 'feature_names_in' in scaler_dict:
            scaler.feature_names_in_ = np.array(scaler_dict['feature_names_in'])
        if 'n_samples_seen' in scaler_dict:
            scaler.n_samples_seen_ = scaler_dict['n_samples_seen']

        artifacts['transforms'][modality] = scaler

    # Reconstruct LabelEncoder/OrdinalEncoder objects
    encoders_json = json_artifacts.get('label_encoders', {})
    for variable, encoder_dict in encoders_json.items():
        if encoder_dict is None:
            artifacts['label_encoders'][variable] = None
            continue

        encoder_type = encoder_dict.get('type')
        
        # Reconstruct LabelEncoder
        if encoder_type == 'LabelEncoder':
            encoder = LabelEncoder()
            encoder.classes_ = np.array(encoder_dict['classes'])
        
        # Reconstruct OrdinalEncoder
        elif encoder_type == 'OrdinalEncoder':
            encoder = OrdinalEncoder(
                handle_unknown=encoder_dict.get('handle_unknown', 'error'),
                unknown_value=encoder_dict.get('unknown_value', None),
                encoded_missing_value=encoder_dict.get('encoded_missing_value', np.nan)
            )
            # Restore fitted attributes
            encoder.categories_ = [np.array(cat) for cat in encoder_dict['categories']]
            if 'n_features_in' in encoder_dict:
                encoder.n_features_in_ = encoder_dict['n_features_in']
            if 'feature_names_in' in encoder_dict:
                encoder.feature_names_in_ = np.array(encoder_dict['feature_names_in'])
            if '_missing_indices' in encoder_dict:
                encoder._missing_indices = encoder_dict['_missing_indices']
            if '_infrequent_enabled' in encoder_dict:
                encoder._infrequent_enabled = encoder_dict['_infrequent_enabled']
        
        else:
            raise ValueError(f"Unknown encoder type: {encoder_type}")

        artifacts['label_encoders'][variable] = encoder

    print(f"[INFO] Loaded artifacts from JSON: {json_path}")
    return artifacts


def convert_joblib_to_json(joblib_path, json_path):
    """
    Convenience function: Convert existing joblib artifacts to JSON.

    Args:
        joblib_path (str): Path to .joblib artifacts file
        json_path (str): Path to save .json artifacts file

    Example:
        >>> convert_joblib_to_json('model.artifacts.joblib', 'model.artifacts.json')
    """
    import joblib
    artifacts = joblib.load(joblib_path)
    artifacts_to_json(artifacts, json_path)
    print(f"[INFO] Converted {joblib_path} → {json_path}")


def convert_json_to_joblib(json_path, joblib_path):
    """
    Convenience function: Convert JSON artifacts back to joblib format.

    Args:
        json_path (str): Path to .json artifacts file
        joblib_path (str): Path to save .joblib artifacts file

    Example:
        >>> convert_json_to_joblib('model.artifacts.json', 'model.artifacts.joblib')
    """
    import joblib
    artifacts = json_to_artifacts(json_path)
    joblib.dump(artifacts, joblib_path)
    print(f"[INFO] Converted {json_path} → {joblib_path}")


if __name__ == '__main__':
    import argparse
    import os

    parser = argparse.ArgumentParser(
        description='Convert flexynesis artifacts between joblib and JSON formats'
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Input file (.joblib or .json)')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Output file (.json or .joblib)')

    args = parser.parse_args()

    # Detect conversion direction based on file extensions
    input_ext = os.path.splitext(args.input)[1].lower()
    output_ext = os.path.splitext(args.output)[1].lower()

    if input_ext in ['.joblib', '.pkl'] and output_ext == '.json':
        # joblib → JSON
        convert_joblib_to_json(args.input, args.output)
    elif input_ext == '.json' and output_ext in ['.joblib', '.pkl']:
        # JSON → joblib
        convert_json_to_joblib(args.input, args.output)
    else:
        print(f"[ERROR] Invalid conversion: {input_ext} → {output_ext}")
        print("[INFO] Supported conversions:")
        print("  - .joblib/.pkl → .json")
        print("  - .json → .joblib/.pkl")
        parser.print_help()
