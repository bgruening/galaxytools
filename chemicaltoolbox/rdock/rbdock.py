import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser(description='Simple wrapper for rbdock')
    parser.add_argument('-n', '--num', type=int, help='Equilibrated system as input')
    parser.add_argument('-s', '--seed', type=int, help='Random seed')
    args = parser.parse_args()

    cmd = ['rbdock', '-i', 'ligands.sdf', '-r', 'receptor.prm', '-p', 'dock.prm', '-n', str(args.num), '-o', 'rdock_output']
    if args.seed != None:
        cmd += ['-s', str(args.seed)]

    ps = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    error_counter = 0
    for stdout_line in iter(ps.stdout.readline, ''):
        print(stdout_line)
        if 'RBT_DOCKING_ERROR' in str(stdout_line):
                error_counter += 1
                if error_counter == 10:
                    print(ps.stdout)
                    exit(23)
        if ps.poll() != None:
            print(ps.stdout)
            exit(int(ps.poll()))

if __name__ == "__main__":
    main()
