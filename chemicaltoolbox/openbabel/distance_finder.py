# Reports distances of ligands to reference points. An example input for the points is:
#
# 5.655   1.497  18.223
# 1.494  -8.367  18.574
# 13.034   6.306  25.232
#
# Data can be space or tab separated but must contain 3 and only 3 numbers for the x, y and z coordinates
#
# That would encode 3 points.
# Each record in the SDF input is read and the closest heavy atom to each of the reference points is recorded as
# a property named distance1 where the numeric part is the index (starting from 1) of the points (in that example
# there would be properties for distance1, distance2 and distance3.

import argparse, os, sys, math
from openbabel import pybel



def log(*args, **kwargs):
    """Log output to STDERR
    """
    print(*args, file=sys.stderr, ** kwargs)


def execute(ligands_sdf, points_file, outfile):
    """
    :param ligands_sdf: A SDF with the 3D molecules to test
    :param points_file: A file with the points to consider.
    :param outfile: The name of the file for the SDF output
    :return:
    """


    points = []

    # read the points
    with open(points_file, 'r') as f:
        for line in f.readlines():
            line.strip()
            if line:
                p = line.split()
                if len(p) == 3:
                    points.append((float(p[0]), float(p[1]), float(p[2])))
                    log("Read points",p)
                    continue
            log("Failed to read line:", line)
    log('Found', len(points), 'atom points')

    sdf_writer = pybel.Outputfile("sdf", outfile, overwrite=True)

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

            p = 0
            for point in points:
                p += 1
                distances = []
                for i in coords:
                    # calculates distance based on cartesian coordinates
                    distance = math.sqrt((point[0] - i[0])**2 + (point[1] - i[1])**2 + (point[2] - i[2])**2)
                    distances.append(distance)
                    # log("distance:", distance)
                min_distance = min(distances)
                # log('Min:', min_distance)
                # log(count, p, min_distance)

                mol.data['distance' + str(p)] = min_distance

            sdf_writer.write(mol)

        except Exception as e:
            log('Failed to handle molecule: '+ str(e))
            continue

    sdf_writer.close()
    log('Wrote', count, 'molecules')


def main():
    global work_dir

    parser = argparse.ArgumentParser(description='XChem distances - measure distances to particular points')

    parser.add_argument('-i', '--input', help="SDF containing the 3D molecules to score)")
    parser.add_argument('-p', '--points', help="PDB format file with atoms")
    parser.add_argument('-o', '--outfile', default='output.sdf', help="File name for results")


    args = parser.parse_args()
    log("XChem distances args: ", args)

    execute(args.input, args.points, args.outfile)


if __name__ == "__main__":
    main()
