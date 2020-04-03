import subprocess
import sys
import os

try:
    int(sys.argv[1])
    num = sys.argv[1]
except (ValueError, IndexError) as e:
    exit(1)

ps = subprocess.Popen(['rbdock', '-i', 'ligands.sdf', '-r', 'receptor.prm', '-p', 'dock.prm', '-n', str(num), '-o', 'output'], stdout=subprocess.PIPE)

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