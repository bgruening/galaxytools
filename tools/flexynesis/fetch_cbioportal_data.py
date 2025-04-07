#!/usr/bin/env python

import argparse
import os

from flexynesis.utils import CBioPortalData


def main():
    parser = argparse.ArgumentParser(description="Fetch and prepare cBioPortal data for Flexynesis.")
    parser.add_argument("--study_id", required=True, help="cBioPortal study ID (e.g., 'brca_tcga')")
    parser.add_argument("--data_types", required=True, help="Comma-separated list of data types (e.g., 'clin,mut,omics')")
    parser.add_argument("--mapped_files", default=None, help="Comma-separated list of .txt files to map to data_types (optional)")
    parser.add_argument("--split_ratio", type=float, default=0.7, help="Training/test split ratio (0.0 to 1.0)")
    parser.add_argument("--output_dir", required=True, help="Output directory for datasets")
    
    args = parser.parse_args()
    
    data_types = args.data_types.split(",")
    if "clin" not in data_types:
        raise ValueError("Clinical data ('clin') is required for splitting the dataset.")
    
    file_mapping = {
        "clin": "data_clinical_patient.txt",  # can be any with 'clinical' in file name
        "mut": "data_mutations.txt",  # any with 'mutations' in file name
        "omics": "data_cna.txt",
        "other": None
    }

    if args.mapped_files:
        mapped_files = args.mapped_files.split(",")
        if len(mapped_files) != len(data_types):
            raise ValueError(f"Number of mapped files ({len(mapped_files)}) must match number of data types ({len(data_types)}).")
        files_to_fetch = {dt: mf for dt, mf in zip(data_types, mapped_files)}
        for mf in mapped_files:
            if not mf.endswith(".txt"):
                raise ValueError(f"Mapped file '{mf}' must end with '.txt'.")
    else:
        files_to_fetch = {dt: file_mapping[dt] for dt in data_types if dt in file_mapping}
    
    invalid_types = set(data_types) - set(file_mapping.keys())
    if invalid_types:
        raise ValueError(f"Invalid data types: {invalid_types}. Supported types: {list(file_mapping.keys())}")

    cbioportal = CBioPortalData(study_id=args.study_id)
    cbioportal.get_cbioportal_data(study_id=args.study_id, files=files_to_fetch)
    dataset = cbioportal.split_data(ratio=args.split_ratio)

    os.makedirs(args.output_dir, exist_ok=True)

    for data_type in data_types:
        if data_type in dataset['train']:
            train_file = os.path.join(args.output_dir, f"{data_type}_train.csv")
            dataset['train'][data_type].to_csv(train_file, index=True)
        if data_type in dataset['test']:
            test_file = os.path.join(args.output_dir, f"{data_type}_test.csv")
            dataset['test'][data_type].to_csv(test_file, index=True)


if __name__ == "__main__":
    main()
