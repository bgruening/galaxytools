# This is a wrapper for JAMM to use it in Galaxy
# It will be accompanied with a jamm.xml file, which specifies the interface that Galaxy is showing to the user
# as well as the JAMM software including JAMM.sh, the bash script that is actually called.

# This wrapper does the following things:
# map the files (from Galaxy history) to a directory
# pass the parameters from the GUI to JAMM.sh
# call JAMM.sh
# map the resulting tabular files back to history items

#import optparse
import argparse, os, shutil, subprocess, sys, tempfile
import shlex
# importing some of the modules used, especially for symlinking, dir manipulation and tool calling.
# since python2.7, argparse is preferred over optparse.

# for reference wrappers, see
# https://bitbucket.org/fubar/rossdev/src/4a91f99b5e1953270c9b35d2ca70c325a10fcfcb/rgdev/bwa_wrapper.py?at=default
# https://www.e-biogenouest.org/wiki/WrapperpythonforGalaxy:Asyntaxexample

def main(): 

    #Read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest = 'input', nargs='+') # + allows for 1 or more arguments
    # parser.add_argument('-c', dest = 'control', nargs='*', default=None) # allow for not supplying -c
    parser.add_argument('-c', dest = 'control', action='append') # allow for not supplying -c
    parser.add_argument('-g', dest = 'gsize')
    parser.add_argument('-o', dest = 'peakfile')
    parser.add_argument('-of', dest = 'filteredpeakfile')
    
    parser.add_argument('-m', dest = 'mode')
    parser.add_argument('-r', dest = 'resolution')
    parser.add_argument('-p', dest = 'processes')
    parser.add_argument('-t', dest = 'type')
    parser.add_argument('-f', dest = 'fraglen')
    parser.add_argument('-b', dest = 'binsize')
    
    args = parser.parse_args()
   
    print "################################################" 
    print "Wrapper debugging" 
    print "################################################" 
    print "Sample files:"
    for j in args.input:
        print j
    if args.control is not None:
        print "Control files:"
        for j in args.control:
            print j

    print "output files:"
    print args.peakfile
    print args.filteredpeakfile
    print "current working dir:"
    print os.getcwd() 
    print "dir with jammwrapper in it:"
    path = os.path.abspath(os.path.dirname(sys.argv[0]))
    print path

    # optparse was depracted, can still be found in may of the example wrappers
    #parser = optparse.OptionParser()
    #parser.add_option( '-i',  dest='input', help='input bed files' )
    #parser.add_option( '-c',  dest='csize', help='chr sizes' )
    #(options, args) = parser.parse_args()	
    
    # create temp dir
    tmp_dir = tempfile.mkdtemp()
    os.mkdir(tmp_dir + "/sample")
    # symlink creation
    for file in args.input:
        filen =  tmp_dir + "/sample/" + os.path.basename(os.path.splitext(file)[0])+".bed"
        print "input files mapped: %s" % filen
        os.symlink(file, filen)
   
    # Here comes some unnecessary repetition
    # if control files are supplied, we make another tmp subdir and put those files there.
    # JAMM.sh is then called with one additional switch -c
    if args.control is not None:
        # symlink creation
        os.mkdir(tmp_dir + "/control")
        for file in args.control:
            filen =  tmp_dir + "/control/" + os.path.basename(os.path.splitext(file)[0])+".bed"
            print "input files mapped: %s" % filen
            os.symlink(file, filen)
        command = ( "bash %s/JAMM.sh -s %s/sample -c %s/control -g %s -o results -m %s -r %s -p %s -t %s -f %s -b %s"
         % ( path, tmp_dir, tmp_dir, args.gsize, args.mode, args.resolution, args.processes, \
             args.type, args.fraglen, args.binsize ) ) 
    else:
        command = ( "bash %s/JAMM.sh -s %s/sample -g %s -o results -m %s -r %s -p %s -t %s -f %s -b %s"
         % ( path, tmp_dir, args.gsize, args.mode, args.resolution, args.processes, \
             args.type, args.fraglen, args.binsize ) ) 

         
    print "Command called by bash:"     
    print command
    # depending on how your programm is called, it may be necessary to use shlex.split on the command string before
    # in this case, this was actually harmful. idk why
#    command = shlex.split(command)
    # Please read the docs for subprocess.Popen. You can't just pass as string as variables need to be extended
    # Note that shell = True can be a security hazard. 
#    proc = subprocess.Popen( command, executable = '/bin/bash', shell=True )
    proc = subprocess.Popen( command, shell=True)
    returncode = proc.wait()
    
    #mv files to a place where galaxy wrapper can find them
    # mvcommand = "mv %s/results/peaks/all.peaks.narrowPeak %s" % ( tmp_dir, args.peakfile )
    mvcommand = "mv results/peaks/all.peaks.narrowPeak %s" % args.peakfile
    os.system(mvcommand)
    mvcommand = "mv results/peaks/filtered.peaks.narrowPeak %s" % args.filteredpeakfile
    os.system(mvcommand)

# clean up temp dir
    if os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )    
    
if __name__ == "__main__":
    main()

