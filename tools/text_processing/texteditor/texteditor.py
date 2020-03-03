import argparse


def filename_to_str(fname, slice=None):
    with open(fname) as f:
        if not slice:
            return f.read()
        return '\n'.join(f.read().split('\n')[slice])

def insert_str_at_pos(str1, str2, pos, line=None):
    if not line:
        return str1[:pos] + str2 + str1[pos:]
    k = str1.split('\n')
    k[line] = insert_str_at_pos(k[line], str2, pos)
    return '\n'.join(k)

def replace_str(str1, substr1, str2):
    return str2.join(str1.split(substr1))

def replace_pos(str1, str2, slice):
    k = str1.split('\n')
    return '\n'.join(k[:slice[0]] + [str2] +  k[slice[1]:])

def slice_file(str, line1, line2):
    if line2 == 0:
        line2 = None
    k = str.split('\n')
    return '\n'.join(k[line1:line2]) + '\n'

def write_output(string, fname):
    with open(fname, 'w') as f:
        f.write(string)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--job', help='insert, replace or delete')
    parser.add_argument('--file1', help='Main file')
    parser.add_argument('--file2', help='File to insert', )
    parser.add_argument('--file2_line1', help='First line of file to insert', type=int)
    parser.add_argument('--file2_line2', help='Last line of file to insert', type=int)
    parser.add_argument('--string_i', help='String to insert', )
    parser.add_argument('--string_d', help='String to delete', )
    parser.add_argument('--pos1', help='First position', type=int)
    parser.add_argument('--line1', help='First line', type=int)
    parser.add_argument('--line2', help='Second line', type=int)
    parser.add_argument('--output', help='Output')

    args = parser.parse_args()


    string = filename_to_str(args.file1)

    if args.file2:
        args.string_i = slice_file(filename_to_str(args.file2), args.file2_line1, args.file2_line2)

    if args.job == 'replace':
        if args.string_d:
            write_output(replace_str(string, args.string_d, args.string_i), args.output)
        elif args.line1 != None and args.line2 != None:
            write_output(replace_pos(string, args.string_i, (args.line1, args.line2)), args.output)
        else:
            raise IOError
    
    elif args.job == 'insert':
        write_output(insert_str_at_pos(string, args.string_i, args.pos1, args.line1), args.output)

    elif args.job == 'delete':
        # just the same as replace but with string_i as an empty string
        if args.string_d:
            write_output(replace_str(string, args.string_d, ''), args.output)
        elif args.line1 != None and args.line2 != None:
            write_output(replace_pos(string, '', (args.line1, args.line2)), args.output)
        else:
            raise IOError
    
    else:
        raise IOError

    # if args.string_i:
        # insert or replace
        
            

if __name__ == "__main__":
    main()


# python texteditor.py --job replace --file1 test-data/t.txt --file2 test-data/s.txt --string_d 'test2' --output test-data/out1.txt
# python texteditor.py --job insert --file1 test-data/t.txt --string_i ins --line1 3 --pos1 2 --output test-data/out2.txt
# python texteditor.py --job insert --file1 test-data/t.txt --string_i ins --line1 3 --pos1 2 --output test-data/out3.txt
# python texteditor.py --job delete --file1 test-data/t.txt --string_d test3  --output test-data/out4.txt
# python texteditor.py --job delete --file1 test-data/t.txt --line1 3 --line2 6  --output test-data/out5.txt