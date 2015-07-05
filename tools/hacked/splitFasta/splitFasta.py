import sys, os
if __name__ == "__main__":
    #assuming perfect input, read every two lines
    inpath = sys.argv[1]
    file_contents = open(inpath, 'r').readlines()
    os.makedirs('splits')
    inname = os.path.basename(inpath)
    for i in range(0, len(file_contents), 2):
        headline = file_contents[i]
        outname = headline[1:headline.index(' ')]+'.fa'
        outfile = open(os.path.join('splits',outname), 'w')
        outfile.write(file_contents[i])
        outfile.write(file_contents[i+1])
        outfile.close()
