import argparse
import json
import numpy as np

def validate_percentages(train_percent, test_percent=None, val_percent=None):
    if test_percent is not None and val_percent is not None:
        total = train_percent + test_percent + val_percent
        if not np.isclose(total, 1.0, atol=1e-5):
            raise ValueError(f"Percentages must sum to 1. Current sum: {total}")
    elif train_percent >= 1:
        raise ValueError(f"Training percentage must be less than 1. Current value: {train_percent}")

def load_data(file_path, file_format):
    print(f"Loading data from {file_path} in {file_format} format")
    # For testing, we'll just return a dummy dataset
    return np.random.rand(1000, 10)

def split_data(data, train_percent, test_percent=None, val_percent=None):
    if test_percent is None:
        # Two file mode
        split_index = int(len(data) * train_percent)
        return data[:split_index], data[split_index:]
    else:
        # Single file mode
        train_end = int(len(data) * train_percent)
        test_end = train_end + int(len(data) * test_percent)
        return data[:train_end], data[train_end:test_end], data[test_end:]

def parse_optuna_param(param_str):
    try:
        return json.loads(param_str)
    except json.JSONDecodeError:
        return param_str

def process_manual_hyperparameters(args):
    return {
        'lr': float(args.lr),
        'batch': int(args.batch),
        'accum': int(args.accum),
        'dropout': float(args.dropout),
        'weight_decay': float(args.weight_decay),
        'warmup_pct': float(args.warmup_pct),
        'lora_rank': int(args.lora_rank),
        'lora_init_scale': float(args.lora_init_scale),
        'lora_scaling_rank': int(args.lora_scaling_rank),
    }

def process_optuna_hyperparameters(args):
    return {
        'lr': parse_optuna_param(args.lr),
        'batch': parse_optuna_param(args.batch),
        'accum': parse_optuna_param(args.accum),
        'dropout': parse_optuna_param(args.dropout),
        'weight_decay': parse_optuna_param(args.weight_decay),
        'warmup_pct': parse_optuna_param(args.warmup_pct),
        'lora_rank': parse_optuna_param(args.lora_rank),
        'lora_init_scale': parse_optuna_param(args.lora_init_scale),
        'lora_scaling_rank': parse_optuna_param(args.lora_scaling_rank),
    }

def main():
    parser = argparse.ArgumentParser(description='Process and validate input data and hyperparameters')
    parser.add_argument('--format', choices=['csv', 'fasta'], required=True)
    parser.add_argument('--input', help='Input file for single file mode')
    parser.add_argument('--train_val', help='Train/val file for two file mode')
    parser.add_argument('--test', help='Test file for two file mode')
    parser.add_argument('--train_percent', type=float, required=True)
    parser.add_argument('--test_percent', type=float)
    parser.add_argument('--val_percent', type=float)
    parser.add_argument('--percentage_summary', required=True, help='Output file for percentage summary')
    
    parser.add_argument('--hyperparameter_method', choices=['manual', 'optuna'], required=True)
    parser.add_argument('--n_trials', type=int, default=10)
    parser.add_argument('--lr', type=str, required=True)
    parser.add_argument('--batch', type=str, required=True)
    parser.add_argument('--accum', type=str, required=True)
    parser.add_argument('--dropout', type=str, required=True)
    parser.add_argument('--weight_decay', type=str, required=True)
    parser.add_argument('--warmup_pct', type=str, required=True)
    parser.add_argument('--lora_rank', type=str, required=True)
    parser.add_argument('--lora_init_scale', type=str, required=True)
    parser.add_argument('--lora_scaling_rank', type=str, required=True)

    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')

    args = parser.parse_args()

    try:
        # Validate percentages
        validate_percentages(args.train_percent, args.test_percent, args.val_percent)

        # Load and split data
        if args.input:
            data = load_data(args.input, args.format)
            train_data, test_data, val_data = split_data(data, args.train_percent, args.test_percent, args.val_percent)
            print(f"Single file mode: Train size: {len(train_data)}, Test size: {len(test_data)}, Val size: {len(val_data)}")
        else:
            train_val_data = load_data(args.train_val, args.format)
            test_data = load_data(args.test, args.format)
            train_data, val_data = split_data(train_val_data, args.train_percent)
            print(f"Two file mode: Train size: {len(train_data)}, Val size: {len(val_data)}, Test size: {len(test_data)}")

        # Process hyperparameters
        if args.hyperparameter_method == 'manual':
            hyperparameters = process_manual_hyperparameters(args)
            print("Manual hyperparameters:")
        else:
            hyperparameters = process_optuna_hyperparameters(args)
            print("Optuna hyperparameter ranges:")
        
        for key, value in hyperparameters.items():
            print(f"  {key}: {value}")

        # Write percentage summary
        with open(args.percentage_summary, 'w') as f:
            if args.input:
                f.write(f"Single file mode:\nTrain: {args.train_percent}\nTest: {args.test_percent}\nValidation: {args.val_percent}\n")
            else:
                f.write(f"Two file mode:\nTrain: {args.train_percent}\nValidation: {1 - args.train_percent}\n")

        print("Data processing and hyperparameter parsing completed successfully.")

    except ValueError as e:
        print(f"Error: {str(e)}")
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")

if __name__ == "__main__":
    main()