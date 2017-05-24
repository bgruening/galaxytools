import os
import re
import sys

shape_file = sys.argv[1]
win_size = int(sys.argv[2])

pattern = re.compile("^>.*$")
toWrite = ""

count_for_id = 1
seq_counter = 0
new_id = ""

seq_id = []
seq_string = []
orig_id = []
name_file = "FASTA/data.names"
array_all_chunks = []
with open(name_file, 'r') as f:
    content = f.read()
    lines = content.split('\n')[:-1]
    for line in lines:
        seq_id.append(int(line.split()[0]))
        seq_string.append(line.split()[1])
        orig_id_srt = line.split()[3]
        orig_id_srt = orig_id_srt.rsplit('_',1)[0]
        orig_id.append(orig_id_srt)


react_dict = {}
react_arr = []

with open(shape_file, 'r') as shape:
    content = shape.read()
    lines = content.split('\n')
    for line in lines:
        if pattern.match(line):
            line = line.replace('>','').strip()
            react_arr=[]
            react_dict[line] = react_arr
            continue
        else:
            react_arr.append(line)

toWrite = ""
chunks = []
for i in range(len(orig_id)):
    if not orig_id[i] in react_dict:
        raise RuntimeError('Error key {} not found'.format(orig_id))

    react_val = react_dict[orig_id[i]]
    toWrite += '>' + str(seq_id[i]) + " " + seq_string[i] + "\n"
    chunks = re.findall(r'\d+', seq_string[i])

    for j in react_val[int(chunks[1])-1:int(chunks[2])]:
        id_s = int(j.split()[0])


        id_s = id_s - int(chunks[1]) + 1
        toWrite += str(id_s) + '\t' + j.split()[1] + "\n"

with open("shape_data_split.react", 'w') as out:
    out.write(toWrite)
