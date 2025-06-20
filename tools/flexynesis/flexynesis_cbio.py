#!/usr/bin/env python

import argparse
import sys

from flexynesis.utils import CBioPortalData


def main():
    parser = argparse.ArgumentParser(description="Fetch and prepare cBioPortal data for Flexynesis.")
    parser.add_argument("--study_id", required=True, help="cBioPortal study ID (e.g., 'brca_tcga')")
    parser.add_argument("--data_types", required=True, help="Comma-separated list of data types (e.g., 'clin,mut,exp')")
    parser.add_argument("--mapped_files", required=True, help="Comma-separated list of .txt files to map to data_types (e.g., 'data_clinical_patient.txt,data_mutations.txt,data_mrna_seq_v2_rsem.txt')")
    parser.add_argument("--split_ratio", type=float, default=0.7, help="Training/test split ratio (0.0 to 1.0)")
    parser.add_argument("--output_dir", required=True, help="Output directory for datasets")

    args = parser.parse_args()

    # Parse and validate inputs
    data_types = [dt.strip() for dt in args.data_types.split(",")]
    mapped_files = [mf.strip() for mf in args.mapped_files.split(",")]

    # Validate that clinical data is included
    if "clin" not in data_types:
        raise ValueError("Clinical data ('clin') is required for splitting the dataset.")

    # Validate that number of data types matches number of files
    if len(data_types) != len(mapped_files):
        raise ValueError(f"Number of data types ({len(data_types)}) must match number of mapped files ({len(mapped_files)}).")

    # Create the files dictionary by zipping data_types and mapped_files
    files = {dt: mf for dt, mf in zip(data_types, mapped_files)}

    # Print configuration for user verification
    print(f"Study ID: {args.study_id}")
    print(f"Split ratio: {args.split_ratio}")
    print("Files mapping:")
    for dt, mf in files.items():
        print(f"  {dt}: {mf}")
    print()

    # Initialize CBioPortalData and fetch data
    try:
        cbioportal = CBioPortalData(study_id=args.study_id)
        print("Downloading and processing data...")
        cbioportal.get_cbioportal_data(study_id=args.study_id, files=files)

        # Check if data was loaded successfully
        if cbioportal.data is None:
            raise ValueError("Error: No data was loaded. Please check the study ID and file mappings.")

        print("Data loaded successfully. Shapes:")
        for key, df in cbioportal.data.items():
            print(f"  {key}: {df.shape}")
        print()

        # Split the data
        print(f"Splitting data with ratio {args.split_ratio}...")
        dataset = cbioportal.split_data(ratio=args.split_ratio)

        # Print split information
        print("Split completed. Dataset shapes:")
        for split in ['train', 'test']:
            print(f"  {split.upper()}:")
            for key, df in dataset[split].items():
                print(f"    {key}: {df.shape}")
        print()

        # Save the datasets
        cbioportal.print_dataset(dataset, outdir=args.output_dir)
        print("Data fetching and preparation completed successfully!")

    except Exception as e:
        print(f"Error occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
