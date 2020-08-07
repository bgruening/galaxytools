from rdkit import Chem
from copy import copy

from pipelines_utils import utils

from molvs import enumerate_tautomers_smiles, canonicalize_tautomer_smiles, Standardizer
from molvs.charge import Uncharger, Reionizer
from standardiser import standardise

standardizer = Standardizer()

def _spam(n):
    out=[]
    for perm in _getPerms(n):
        elem = [ int(i) for i in list(perm) ]
        out.append(elem)
    return out

def _getPerms(n):
    from itertools import permutations
    for i in _getCandidates(n):
        for perm in set(permutations(i)):
            yield ''.join(perm)

def _getCandidates(n):
    for i in range(0, n+1):
        res = "1" * i + "0" * (n - i)
        yield res

def enumerateTautomers(mol):
    """
    Get all of the Tautomers of a given molecule
    :param mol: the input molecule
    :return: a list of Tautomers
    """
    smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
    tauts = enumerate_tautomers_smiles(smiles)
    ##TODO Append Parent molecule name
    return  [Chem.MolFromSmiles(x) for x in tauts]

def getCanonTautomer(mol):
    """
    Get the canonical tautomer form
    :param mol: the input molecule
    :return: a list of Tautomers
    """
    smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
    x = canonicalize_tautomer_smiles(smiles)
    return Chem.MolFromSmiles(x)


def enumerateStereoIsomers(mol):
    out = []
    chiralCentres = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    #return the molecule object when no chiral centres where identified
    if chiralCentres == []:
        return [mol]

    #All bit permutations with number of bits equals number of chiralCentres
    elements = _spam(len(chiralCentres))

    for isoId,element in enumerate(elements):
        for centreId,i in enumerate(element):
            atomId = chiralCentres[centreId][0]
            if i == 0:
                mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
            elif i == 1:
                mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        outmol = copy(mol)
        utils.log("Enumerated ", Chem.MolToSmiles(mol, isomericSmiles=True))
        out.append(outmol)
    return out


def molVsStandardizer(mol):
    return standardizer.standardize(mol)

def flatkinsonStandardizer(mol):
    return standardise.run(mol)

STANDARD_MOL_METHODS = {"molvs": molVsStandardizer, "flatkinson": flatkinsonStandardizer}

def getNeutralMolecule(mol):
    uncharger = Uncharger()
    return uncharger.uncharge(mol)

def getReionisedMolecule(mol):
    reioniser = Reionizer()
    return reioniser.reionize(mol)