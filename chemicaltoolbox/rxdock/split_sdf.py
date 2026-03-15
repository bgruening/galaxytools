import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Chunk SDF for parallelizing rDock.")
    parser.add_argument(
        "-c", "--chunks", type=int, help="Number of chunks to generate"
    )
    parser.add_argument("-n", "--number", type=int, help="Number of molecules in the SDF")
    args = parser.parse_args()

    file_size = int(args.number / args.chunks + 1)
    file_index = 0
    molecules_added = 0
    f = open(f"ligands{file_index:04}.sdf", 'wb')

    for line in sys.stdin.buffer.raw:
        f.write(line)
        if b'$$$$' in line:
            molecules_added += 1
            if molecules_added == file_size:
                f.close()
                file_index += 1
                molecules_added = 0
                f = open(f'ligands{file_index:04}.sdf', 'wb')
    f.close()
    if file_index > args.chunks:  # should be impossible
        raise Exception("Too many chunks were created!")


if __name__ == "__main__":
    main()
