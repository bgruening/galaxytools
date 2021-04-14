import argparse

from chembl_structure_pipeline import checker, standardizer


def load_mols(input_file):
    """
    Returns a list of strings, each a molblock
    """
    with open(input_file) as f:
        mols = [
            "".join(("\n", mol.strip())) for mol in f.read().strip().split("$$$$\n")
        ]
    return mols


def write_mols(mols, output_file):
    """
    Writes a list of molblocks to an SDF
    """
    with open(output_file, "w") as f:
        f.write("\n$$$$".join(mols))


def standardize_molblock(mol):
    return standardizer.standardize_molblock(mol)


def get_parent_molblock(mol):
    return standardizer.get_parent_molblock(mol)[0]


def check_molblock(mol):
    issues = checker.check_molblock(mol)
    max_penalty_score = str(max([issue[0] for issue in issues])) if issues else "0"
    message = "; ".join([issue[1] for issue in issues])
    mol_with_issues = "\n".join(
        (mol, "> <MaxPenaltyScore>", max_penalty_score, "> <IssueMessages>", message)
    )
    return mol_with_issues


def main():
    parser = argparse.ArgumentParser(description="Search ChEMBL database for compounds")
    parser.add_argument("-i", "--input", help="SDF/MOL input")
    parser.add_argument("-o", "--output", help="Standardized output")
    parser.add_argument(
        "--standardize", action="store_true", help="Standardize molblock"
    )
    parser.add_argument(
        "--get_parent", action="store_true", help="Get parent molblock."
    )
    parser.add_argument("--check", action="store_true", help="Check molblock")
    args = parser.parse_args()

    mols = load_mols(args.input)
    if args.standardize:
        mols = [standardize_molblock(mol) for mol in mols]
    if args.get_parent:
        mols = [get_parent_molblock(mol) for mol in mols]
    if args.check:
        mols = [check_molblock(mol) for mol in mols]
    write_mols(mols, args.output)


if __name__ == "__main__":
    main()
