# Create dir containing ligands.sdf and protein.pdb
# Enter docker container like this:
#   docker run -it --rm --gpus all -v $PWD:/root/train/fragalysis_test_files/work:Z informaticsmatters/deep-app-ubuntu-1604:latest bash
#
# Now inside the container run like this:
#   mkdir /tmp/work
#   rm -rf /tmp/work/* && python3 work/transfs.py -i work/test-data/ligands.sdf -r work/test-data/receptor.pdb -w /tmp/work
#
# If testing with no GPU you can use the --mock option to generate random scores
#
# Start container for testing like this:
#    docker run -it --rm -v $PWD:$PWD:Z -w $PWD informaticsmatters/deep-app-ubuntu-1604:latest bash
# Inside container test like this:
#   mkdir /tmp/work
#   cd chemicaltoolbox/xchem-deep
#   rm -rf /tmp/work/* && python3 transfs.py -i test-data/ligands.sdf -r test-data/receptor.pdb -w /tmp/work --mock
#

import argparse, os, sys, math
import subprocess
import random
from openbabel import pybel

types_file_name = 'inputs.types'
types_file_name = 'inputs.types'
predict_file_name = 'predictions.txt'
work_dir = '.'
paths = None
inputs_protein = []
inputs_ligands = []


def log(*args, **kwargs):
    """Log output to STDERR
    """
    print(*args, file=sys.stderr, ** kwargs)

def write_raw_inputs(receptor_pdb, ligands_sdf, distance):

    global work_dir
    global inputs_protein
    global inputs_ligands
    global paths


    log("Writing data to", work_dir)
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)

    receptor_file = os.path.basename(receptor_pdb)

    sdf_writers = {}
    paths = []

    with open(receptor_pdb, 'r') as f:
        lines = f.readlines()

    count = 0
    for mol in pybel.readfile("sdf", ligands_sdf):
        count += 1
        if count % 50000 == 0:
            log('Processed', count)

        try:
            # print("Processing mol", mol.title)

            clone = pybel.Molecule(mol)
            clone.removeh()

            coords = []
            for atom in clone.atoms:
                coords.append(atom.coords)

            watnumcode = ''

            # getting receptor without waters that will clash with ligand
            new_receptor_pdb = []
            for line in lines:
                if line[17:20] == 'HOH':
                    x, y, z = float(line[30:39]),  float(line[39:46]), float(line[46:55])
                    distances = []
                    for i in coords:
                        distances.append(math.sqrt((x-i[0])**2 + (y-i[1])**2 + (z-i[2])**2))  # calculates distance based on cartesian coordinates
                    if min(distances) > distance: # if all distances are larger than 2.0A, then molecule makes it to new file
                        new_receptor_pdb.append(line)
                    else:
                        watnum = line[23:28].strip()
                        # print("Skipped water " + watnum)
                        watnumcode += '-' + watnum
                if line[17:20] != 'LIG' and line[17:20] != 'HOH':  # ligand lines are also removed
                    new_receptor_pdb.append(line)


            name = receptor_file[0:-4] + watnumcode
            # print('CODE:', name)
            mol.data['TransFSReceptor'] = name


            if watnumcode not in sdf_writers:

                dir = os.path.sep.join([work_dir, name])
                log('WRITING to :', dir)

                os.mkdir(dir)
                paths.append(dir)

                sdf_writers[watnumcode] = pybel.Outputfile("sdf", os.path.sep.join([dir, 'ligands.sdf']))

                # writing into new pdb file
                receptor_writer = open(os.path.sep.join([dir, 'receptor.pdb']), "w+")
                for line in new_receptor_pdb:
                    receptor_writer.write(str(line))
                receptor_writer.close()

            sdf_out = sdf_writers[watnumcode]
            sdf_out.write(mol)

        except Exception as e:
            log('Failed to handle molecule: '+ str(e))
            continue

    for writer in sdf_writers.values():
        writer.close()

    log('Wrote', count, 'molecules and', len(sdf_writers), 'proteins')

def write_inputs(protein_file, ligands_file, distance):
    global types_file_name
    global work_dir
    global inputs_protein
    global inputs_ligands
    global prepared_ligands

    write_raw_inputs(protein_file, ligands_file, distance)

    types_path = os.path.sep.join([work_dir, types_file_name])
    log("Writing types to", types_path)
    with open(types_path, 'w') as types_file:

        for path in paths:

            log("Gninatyping ligands in", path)

            ligands_dir = os.path.sep.join([path, 'ligands'])
            os.mkdir(ligands_dir)

            cmd1 = ['gninatyper', os.path.sep.join([path, 'ligands.sdf']), os.path.sep.join([ligands_dir, 'ligand'])]
            log('CMD:', cmd1)
            exit_code = subprocess.call(cmd1)
            log("Status:", exit_code)
            if exit_code:
                raise Exception("Failed to write ligands")
            ligand_gninatypes = os.listdir(os.path.sep.join([path, 'ligands']))

            log("Gninatyping proteins in", path)

            proteins_dir = os.path.sep.join([path, 'proteins'])
            os.mkdir(proteins_dir)

            cmd2 = ['gninatyper', os.path.sep.join([path, 'receptor.pdb']), os.path.sep.join([proteins_dir, 'protein'])]
            log('CMD:', cmd2)
            exit_code = subprocess.call(cmd2)
            log("Status:", exit_code)
            if exit_code:
                raise Exception("Failed to write proteins")
            protein_gninatypes = os.listdir(os.path.sep.join([path, 'proteins']))

            num_proteins = 0
            num_ligands = 0

            for protein in protein_gninatypes:
                num_proteins += 1
                num_ligands = 0
                inputs_protein.append(protein)
                inputs_protein.append(os.path.sep.join([path, 'proteins', protein]))
                for ligand in ligand_gninatypes:
                    num_ligands += 1
                    log("Handling", protein, ligand)
                    inputs_ligands.append(os.path.sep.join([path, 'ligands', ligand]))
                    line = "0 {0}{3}proteins{3}{1} {0}{3}ligands{3}{2}\n".format(path, protein, ligand, os.path.sep)
                    types_file.write(line)

    return num_proteins, num_ligands


def generate_predictions_filename(work_dir, predict_file_name):
    return "{0}{1}{2}".format(work_dir, os.path.sep, predict_file_name)


def run_predictions():
    global types_file_name
    global predict_file_name
    global work_dir
    # python3 scripts/predict.py -m resources/dense.prototxt -w resources/weights.caffemodel -i work_0/test_set.types >> work_0/caffe_output/predictions.txt
    cmd1 = ['python3', '/train/fragalysis_test_files/scripts/predict.py',
            '-m', '/train/fragalysis_test_files/resources/dense.prototxt',
            '-w', '/train/fragalysis_test_files/resources/weights.caffemodel',
            '-i', os.path.sep.join([work_dir, types_file_name]),
            '-o', os.path.sep.join([work_dir, predict_file_name])]
    log("CMD:", cmd1)
    subprocess.call(cmd1)


def mock_predictions():
    global work_dir
    global predict_file_name

    log("WARNING: generating mock results instead of running on GPU")
    outfile = generate_predictions_filename(work_dir, predict_file_name)
    count = 0
    with open(outfile, 'w') as predictions:
        for path in paths:
            log("Reading", path)
            protein_gninatypes = os.listdir(os.path.sep.join([path, 'proteins']))
            ligand_gninatypes = os.listdir(os.path.sep.join([path, 'ligands']))
            for protein in protein_gninatypes:
                for ligand in ligand_gninatypes:
                    count += 1
                    score = random.random()
                    line = "{0} | 0 {1}{4}proteins{4}{2} {1}{4}ligands{4}{3}\n".format(score, path, protein, ligand,
                                                                                       os.path.sep)
                    # log("Writing", line)
                    predictions.write(line)

    log('Wrote', count, 'mock predictions')


def read_predictions():
    global predict_file_name
    global work_dir
    scores = {}
    with open("{0}{1}{2}".format(work_dir, os.path.sep, predict_file_name), 'r') as input:
        for line in input:
            # log(line)
            tokens = line.split()
            if len(tokens) == 5 and tokens[1] == '|':
                # log(len(tokens), tokens[0], tokens[3], tokens[4])
                record_no = inputs_ligands.index(tokens[4])
                if record_no is not None:
                    # log(record_no, tokens[0])
                    scores[record_no] = tokens[0]
    log("Found", len(scores), "scores")
    return scores


def patch_scores_sdf(outfile, scores):

    counter = 0
    sdf_path = "{0}{1}{2}".format(work_dir, os.path.sep, outfile)
    log("Writing results to {0}".format(sdf_path))
    sdf_file = pybel.Outputfile("sdf", sdf_path)

    for path in paths:
        for mol in pybel.readfile("sdf", os.path.sep.join([path, 'ligands.sdf'])):
            if counter in scores:
                score = scores[counter]
                # og("Score for record {0} is {1}".format(counter, score))
                mol.data['TransFSScore'] = score
                sdf_file.write(mol)
            else:
                log("No score found for record", counter)
            counter += 1
    sdf_file.close()


def execute(ligands_sdf, protein, outfile, distance, mock=False):

    write_inputs(protein, ligands_sdf, distance)
    if mock:
        mock_predictions()
    else:
        run_predictions()
    scores = read_predictions()
    patch_scores_sdf(outfile, scores)


def main():
    global work_dir

    parser = argparse.ArgumentParser(description='XChem deep - pose scoring')

    parser.add_argument('-i', '--input', help="SDF containing the poses to score)")
    parser.add_argument('-r', '--receptor', help="Receptor file for scoring (PDB or Mol2 format)")
    parser.add_argument('-d', '--distance', type=float, default=2.0, help="Cuttoff for removing waters")
    parser.add_argument('-o', '--outfile', default='output.sdf', help="File name for results")
    parser.add_argument('-w', '--work-dir', default=".", help="Working directory")
    parser.add_argument('--mock', action='store_true', help='Generate mock scores rather than run on GPU')

    args = parser.parse_args()
    log("XChem deep args: ", args)

    work_dir = args.work_dir

    execute(args.input, args.receptor, args.outfile, args.distance, mock=args.mock)


if __name__ == "__main__":
    main()
