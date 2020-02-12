# Create dir containing ligands.sdf and protein.pdb
# Enter docker container like this:
#   docker run -it --rm --gpus all -v $PWD:/root/train/fragalysis_test_files/work:Z informaticsmatters/deep-app-ubuntu-1604:latest bash
#
# Now inside the container run like this:
#   rm -rf work/* && python3 -m  pipelines_obabel.gnina.xchem_deep_score -i ligands.sdf -r protein.pdb -w work
#
# If testing with no GPU you can use the --mock option to generate random scores
#
# Start container for testing like this:
#    docker run -it --rm -v $PWD:$PWD:Z -w $PWD informaticsmatters/deep-app-ubuntu-1604:latest bash
# Inside container test like this:
#   mkdir /tmp/work
#   cd chemicaltoolbox/xchem-deep
#   rm -rf /tmp/work/* && python3 xchem_deep_score.py -i test-data/ligands.sdf -r test-data/receptor.pdb -w /tmp/work --mock
#

import argparse, os, re, sys
import random
from openbabel import pybel

types_file_name = 'inputs.types'
types_file_name = 'inputs.types'
predict_file_name = 'predictions.txt'
work_dir = '.'
inputs_protein = []
inputs_ligands = []

def log(*args, **kwargs):
    """Log output to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)

def write_inputs(protein_file, ligands_file):
    global types_file_name
    global work_dir
    global inputs_protein
    global inputs_ligands

    ligands_path = "{0}{1}ligands".format(work_dir, os.path.sep)
    log("Writing ligands to", ligands_path)
    os.mkdir(ligands_path)
    cmd1 = "gninatyper {0} {1}{2}ligands{2}ligand".format(ligands_file, work_dir, os.path.sep)
    log('CMD:', cmd1)
    exit_code = os.system(cmd1)
    log("Status:", exit_code)
    ligand_gninatypes = os.listdir("{0}{1}ligands".format(work_dir, os.path.sep))

    proteins_path = "{0}{1}proteins".format(work_dir, os.path.sep)
    log("Writing proteins to", proteins_path)
    os.mkdir(proteins_path)
    cmd2 = "gninatyper {0} {1}{2}proteins{2}protein".format(protein_file, work_dir, os.path.sep)
    log('CMD:', cmd2)
    exit_code = os.system(cmd2)
    log("Status:", exit_code)
    protein_gninatypes = os.listdir("{0}{1}proteins".format(work_dir, os.path.sep))

    types_path = "{0}{1}{2}".format(work_dir, os.path.sep, types_file_name)
    log("Writing types to", types_path)
    num_proteins = 0
    num_ligands = 0
    with open(types_path, 'w') as types_file:
        for protein in protein_gninatypes:
            num_proteins += 1
            num_ligands = 0
            inputs_protein.append(protein)
            for ligand in ligand_gninatypes:
                num_ligands += 1
                log("Handling", protein, ligand)
                inputs_ligands.append(ligand)
                line = "0 {0}{3}proteins{3}{1} {0}{3}ligands{3}{2}\n".format(work_dir, protein, ligand, os.path.sep)
                types_file.write(line)

    return num_proteins, num_ligands

def generate_predictions_filename(work_dir, predict_file_name):
    return "{0}{1}{2}".format(work_dir, os.path.sep, predict_file_name)

def run_predictions():
    global types_file_name
    global predict_file_name
    global work_dir
    # python3 scripts/predict.py -m resources/dense.prototxt -w resources/weights.caffemodel -i work_0/test_set.types >> work_0/caffe_output/predictions.txt
    cmd1 = "python3 /root/train/fragalysis_test_files/scripts/predict.py -m /root/train/fragalysis_test_files/resources/dense.prototxt" +\
           " -w /root/train/fragalysis_test_files/resources/weights.caffemodel" +\
            " -i {0}/{1} -o {0}/{2}".format(work_dir, types_file_name, predict_file_name)
    log("CMD:", cmd1)
    os.system(cmd1)

def mock_predictions():
    global work_dir
    global predict_file_name
    global inputs_protein
    global inputs_ligands
    log("WARNING: generating mock results instead of running on GPU")
    outfile = generate_predictions_filename(work_dir, predict_file_name)
    with open(outfile, 'w') as predictions:
        for protein in inputs_protein:
            for ligand in inputs_ligands:
                score = random.random()
                line = "{0} | 0 {1}{4}proteins{4}{2} {1}{4}ligands{4}{3}\n".format(score, work_dir, protein, ligand, os.path.sep)
                # log("Writing", line)
                predictions.write(line)


def read_predictions():
    global predict_file_name
    global work_dir
    scores = {}
    with open("{0}{1}{2}".format(work_dir, os.path.sep, predict_file_name), 'r') as input:
        for line in input:
            #log(line)
            tokens = line.split()
            if len(tokens) == 5 and tokens[1] == '|':
                # log(len(tokens), tokens[0], tokens[3], tokens[4])
                record_no = match_ligand(tokens[4])
                if record_no is not None:
                    # log(record_no, tokens[0])
                    scores[record_no] = tokens[0]
    log("Found", len(scores), "scores")
    return scores

def patch_scores_sdf(sdf_in, outfile, scores):

    global work_dir

    counter = 0
    sdf_path = "{0}{1}{2}.sdf".format(work_dir, os.path.sep, outfile)
    log("Writing results to {0}".format(sdf_path))
    sdf_file = pybel.Outputfile("sdf", sdf_path)
    for mol in pybel.readfile("sdf", sdf_in):
        if counter in scores:
            score = scores[counter]
            # og("Score for record {0} is {1}".format(counter, score))
            mol.data['XChemDeepScore'] = score
            sdf_file.write(mol)
        else:
            log("No score found for record", counter)
        counter += 1
    sdf_file.close()


# work/ligands/ligand_9.gninatypes
ligand_patt = re.compile(r'.*ligands/ligand_(\d+)\.gninatypes$')

def match_ligand(s):
    global ligand_patt
    m = ligand_patt.match(s)
    if m:
        i = m.group(1)
        return int(i)
    else:
        return None


def execute(ligands_sdf, protein, outfile, mock=False):

    output_base = "{0}{1}{2}".format(work_dir, os.path.sep, outfile)

    write_inputs(protein, ligands_sdf)
    if mock:
        mock_predictions()
    else:
        run_predictions()
    scores = read_predictions()
    patch_scores_sdf(ligands_sdf, outfile, scores)

def main():

    global work_dir

    parser = argparse.ArgumentParser(description='XChem deep - pose scoring')

    parser.add_argument('-i', '--input', help="SDF containing the poses to score)")
    parser.add_argument('-r', '--receptor', help="Receptor file for scoring (PDB or Mol2 format)")
    parser.add_argument('-o', '--outfile', default='output', help="Base file name for results")
    parser.add_argument('-w', '--work-dir', default=".", help="Working directory")
    parser.add_argument('--mock', action='store_true', help='Generate mock scores rather than run on GPU')

    args = parser.parse_args()
    log("XChem deep args: ", args)

    work_dir = args.work_dir

    execute(args.input, args.receptor, args.outfile, mock=args.mock)

if __name__ == "__main__":
    main()
