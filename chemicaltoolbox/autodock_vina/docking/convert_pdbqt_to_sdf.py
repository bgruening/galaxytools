import pybel, openbabel
import sys


def main():
	if len(sys.argv) == 3:
		process(sys.argv[1], sys.argv[2])
	else:
		print("Usage: convert_pdbqt_to_sdf.py <input-pdbqt-file> <output-sdf-file>")
		exit(1)

def add_property(mol, prop_name, prop_value):
	newData = openbabel.OBPairData() 
	newData.SetAttribute(prop_name)
	newData.SetValue(prop_value) 
	mol.OBMol.CloneData(newData)

def process(input, output):
	docked = pybel.readfile('pdbqt', input)
	sdf = pybel.Outputfile("sdf", output, overwrite=True)
	for mol in docked:
		if mol.OBMol.HasData('REMARK'):
			remark = mol.OBMol.GetData('REMARK').GetValue()
			lines = remark.splitlines()
			tokens = lines[0].split()
		
			# add the score property
			add_property(mol, "SCORE", tokens[2]) 
			# add the first RMSD property
			add_property(mol, "RMSD_LB", tokens[3])
			# add the second RMSD property
			add_property(mol, "RMSD_UB", tokens[4])

		sdf.write(mol)

	sdf.close()

if __name__ == "__main__":
    main()

