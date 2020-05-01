import argparse

def get_coordinates(lines):
    version = lines[3][34:39]
    molecule = []
    if version == 'V2000':
        natom = int(lines[3][:3].strip())
        for i in range(1, natom + 1):
            temp = []
            j = 3 + i
            x = float(lines[j][:10].strip())
            y = float(lines[j][11:20].strip())
            z = float(lines[j][21:30].strip())
            temp.extend([x, y, z])
            molecule.append(temp)
    else:
        read = 0
        for line in lines:
            if "END ATOM" in line:
                read = 0
                break
            if read:
                temp = []
                a = line.split(" ")
                x, y, z = float(a[5]), float(a[6]), float(a[7])
                temp.extend([x, y, z])
                molecule.append(temp)
            if "BEGIN ATOM" in line:
                read = 1
    return molecule


def select_points(all_coordinates):
    tol = 1.5
    select = []

    for molecule in all_coordinates:
        for coordinates in molecule:
            tv = 0
            temp = []
            x, y, z = coordinates
            for record in select:
                xr, yr, zr = record
                if xr-tol < x and x < xr+tol and \
                   yr-tol < y and y < yr+tol and \
                   zr-tol < z and z < zr+tol:
                    tv = 1
                    break
            if tv == 1:
                continue
            temp.extend([x, y, z])
            select.append(temp)
    return select


def sdfout(centers, writer):
    n = len(centers)
    writer.write("Frankenstein_ligand\nGalaxy select_points_sdf tool\n\n")
    writer.write("%3d  0  0  0  0  0  0  0  0  0999 V2000\n" % n)
    for record in centers:
        x, y, z = record
        writer.write("%10.4f%10.4f%10.4f C   0  0  0  0  0  0  0  0  0  0  0  0\n" % (x, y, z))

    writer.write("M  END\n$$$$\n")


def main():
    parser = argparse.ArgumentParser(description='RDKit screen')
    parser.add_argument('-i', '--input',
                        help="Input file")
    parser.add_argument('-o', '--output',
                        help="Base name for output file (no extension).")
    args = parser.parse_args()

    mol_coordinates = []
    all_coordinates = []
    with open(args.input) as file:
        for line in file:
            if line.strip() == '$$$$':
                temp = get_coordinates(mol_coordinates)
                all_coordinates.append(temp)
                mol_coordinates.clear()
            else:
                mol_coordinates.append(line)
    centers = select_points(all_coordinates)
    with open(args.output, 'w+') as writer:
        sdfout(centers, writer)

if __name__ == "__main__":
    main()
