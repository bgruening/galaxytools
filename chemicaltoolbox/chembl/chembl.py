from chembl_webresource_client.new_client import new_client
import argparse

def open_file(filename):
    with open(filename) as f:
        return f.read().split('\n')[0]

def get_smiles(res):
    """
    Get a list of SMILES from function results
    """ 
    smiles = []
    for smi in res: 
        smiles.append(smi['molecule_structures']['canonical_smiles']) 
    return smiles

def sim_search(smiles, tanimoto):
    """
    Return compounds which are within a Tanimoto range of the SMILES input
    """
    similarity = new_client.similarity
    return similarity.filter(smiles=smiles, similarity=tanimoto).only(['molecule_structures'])
    
def substr_search(smiles):
    """
    Return compounds which contain the SMILES substructure input
    """
    substructure = new_client.substructure
    return substructure.filter(smiles=smiles).only(['molecule_structures'])
    
def filter_drugs(mols):
    """
    Return only compounds which are approved drugs
    """
    return mols.filter(max_phase=4)

def filter_biotherapeutic(mols):
    """
    Return only biotherapeutic molecules
    """
    return mols.filter(biotherapeutic__isnull=False)

def filter_nat_prod(mols):
    """
    Return only natural products
    """
    return mols.filter(natural_product=1)

def filter_ro5(mols):
    """
    Return only compounds with no RO5 violations
    """
    return mols.filter(molecule_properties__num_ro5_violations=0)

def main():
    parser = argparse.ArgumentParser(description='Search ChemBL database for compounds')
    parser.add_argument('-i', '--input', help='SMILES input')
    parser.add_argument('-f', '--file', help='SMILES input as file')
    parser.add_argument('-o', '--output', help="SMILES output")
    parser.add_argument('-t', '--tanimoto', type=int, help='Tanimoto similarity score')
    parser.add_argument('-s', '--substructure', action='store_true', help='Substructure search using the SMILES input.')
    parser.add_argument('-d', '--drugs', action='store_true', help='Filter approved drugs')
    parser.add_argument('-b', '--biotherapeutic', action='store_true', help='Filter biotherapeutic molecules')
    parser.add_argument('-n', '--nat-prod', action='store_true', help='Filter natural products')
    parser.add_argument('-r', '--ro5', action='store_true', help='Filter compounds that pass Lipinski RO5')

    args = parser.parse_args()

    if args.file:  # get SMILES from file rather than -i option
        args.input = open_file(args.file)

    if len(args.input) < 5:
        raise IOError('SMILES must be at least 5 characters long.')

    if args.substructure:  # specify search type: substructure or similarity
        mols = substr_search(args.input)
    else:
        mols = sim_search(args.input, args.tanimoto)

    # filter options:
    if args.drugs:
        mols = filter_drugs(mols)

    if args.biotherapeutic:
        mols = filter_biotherapeutic(mols)

    if args.nat_prod:
        mols = filter_nat_prod(mols)

    if args.ro5:
        mols = filter_ro5(mols)

    # get SMILES from search output
    mols = get_smiles(mols)

    # write to file
    with open(args.output, 'w') as f:
        f.write('\n'.join(mols))
    

if __name__ == "__main__":
    main()
