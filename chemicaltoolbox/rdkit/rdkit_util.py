from rdkit import Chem


def get_supplier(infile, format='smiles'):
    """
    Returns a generator over a SMILES or InChI file. Every element is of RDKit
    molecule and has its original string as _Name property.
    """
    with open(infile) as handle:
        for line in handle:
            line = line.strip()
            if format == 'smiles':
                mol = Chem.MolFromSmiles(line, sanitize=True)
            elif format == 'inchi':
                mol = Chem.inchi.MolFromInchi(line, sanitize=True, removeHs=True, logLevel=None, treatWarningAsError=False)
            if mol is None:
                yield False
            else:
                mol.SetProp('_Name', line.split('\t')[0])
                yield mol
