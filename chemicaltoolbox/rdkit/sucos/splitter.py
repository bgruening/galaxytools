import argparse
from rdkit import Chem
import utils


def process(input, basename, size):

    suppl = Chem.ForwardSDMolSupplier(input)
    count = 0
    chunk = 0
    writer = None
    for mol in suppl:
        if count % size == 0:
            if writer != None:
                writer.close()
            chunk += 1
            writer = new_writer(basename, chunk)
        count +=1
        if mol == None:
            utils.log("Skipping record", count)
            continue
        writer.write(mol)
        writer.flush()
    if writer != None:
        writer.close()



def new_writer(basename, count):
    name = basename + "_" + str(count) + ".sdf"
    utils.log("Creating new output", name)
    writer = Chem.SDWriter(name)
    return writer


def main():

    parser = argparse.ArgumentParser(description='Split SDF with RDKit')
    parser.add_argument('-i', '--input', help='Input file in SDF format.')
    parser.add_argument('-o', '--output', default="output", help='Base name for output files in SDF format.')
    parser.add_argument('-c', '--chunk-size', help='Number of records in each file', type=int, default=10)

    args = parser.parse_args()
    utils.log("Splitter Args: ", args)


    process(args.input, args.output, size=args.chunk_size)


if __name__ == "__main__":
    main()