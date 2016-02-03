#!/usr/bin/python2

def get_argument_parser():
    parser = argparse.ArgumentParser(
        description="RCAS provides intuitive reports and publication ready graphics"
        " from input peak intervals in BED foramt,"
        " which are detected in clip-seq data."
        )

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s 0.1')

    parser.add_argument("BED",
                        metavar="BED",
                        nargs="*",
                        help="Target intervals in BED format.")

    parser.add_argument("--RCAS_path", "-r",
                        metavar="path/to/RCAS",
                        required=True,
                        help="Path to RCAS.")

    parser.add_argument("--genome", "-g",
                        metavar="FILE",
                        required=True,
                        help="Reference genome whose version should conform with"
                        " generation of input peak invervals.")

    parser.add_argument("--gff3", "-f",
                        metavar="FILE",
                        required=True,
                        help="Annotation reference in gff3 format"
                        " whose version should conform with"
                        " genome version.")

    parser.add_argument("--species", "-s",
                        default="human",
                        choices=["human", "fly", "worm", "mouse"],
                        help="Required for running GO-term or pathway enrichment."
                        " Default: human")

    parser.add_argument("--cores", "-j",
                        metavar="N",
                        default="1",
                        help="Use at most N cores in parallel (default: 1).")

    parser.add_argument("--forcerun", "-F",
                        action='store_true',
                        help="Force the re-execution of snakemake."
                              " Use this option if you want to have all"
                              " output in your workflow updated.")

    parser.add_argument("--run_motif", "-m",
                        action='store_true',
                        help="Run motif search.")

    parser.add_argument("--run_PATHrich", "-p",
                        action='store_true',
                        help="Run pathway enrichment.")

    parser.add_argument("--run_GOrich", "-t",
                        action='store_true',
                        help="Run GO term enrichment.")

    parser.add_argument("--run_coverage", "-c",
                        action='store_true',
                        help="Run coverage profile.")

    parser.add_argument("--run_all", "-a",
                        action='store_true',
                        help="Run all steps.")

    return parser

def extract_key(filename):
    base = os.path.basename(filename)
    key = base.split(".")[:-1]
    key = ".".join(key)
    return key

def generate_config(args):
    BED_files = args.BED
    infiles = {}

    for BED in BED_files:
        infiles[extract_key(BED)] = BED

    if args.run_motif:
        run_motif = "True"
    else:
        run_motif = "False"

    if args.run_PATHrich:
        run_PATHrich = "True"
    else:
        run_PATHrich = "False"

    if args.run_GOrich:
        run_GOrich = "True"
    else:
        run_GOrich = "False"

    if args.run_coverage:
        run_coverage = "True"
    else:
        run_coverage = "False"

    if args.run_all:
        run_motif = "True"
        run_PATHrich = "True"
        run_GOrich = "True"
        run_coverage = "True"

    config = {
      "RCAS_path": os.path.abspath(args.RCAS_path),

      "gff3": args.gff3,

      "genome": args.genome,

      "species": args.species,

      "infile": infiles,

      "switch": {
        "run_motif": run_motif,
        "run_PATHrich": run_PATHrich,
        "run_GOrich": run_GOrich,
        "run_coverage": run_coverage
      }
    }

    # Writing JSON data
    with open('config.json', 'w') as f:
         json.dump(config, f, sort_keys=True, indent=4)

    print "\nwrote config.json.\n"

def call_snakemake(RCAS_path, forcerun, cores):
    print "start snakemake:\n"

    if forcerun:
        forcerun = "F"
    else:
        forcerun = ""

    command_line = "snakemake -%sp -j %s -s %s/src/RCAS.snakefile" % (forcerun, cores, RCAS_path)
    cmd = shlex.split(command_line)

    p = Popen(cmd)

    try:
        p.communicate()
    except KeyboardInterrupt:
        try:
            p.send_signal(signal.SIGINT)
            time.sleep(0.1)
            if p.poll() is not None:
                print "Process is interrupted."
                return
            p.terminate()
            time.sleep(0.1)
            if p.poll() is not None:
                print "Process is terminated."
                return
            p.kill()
            print "Process is killed."
        except OSError:
            pass
        except Exception as e:
            print "Error while terminating subprocess (pid=%i): %s" \
                % (p.pid, e)
        return

if __name__ == '__main__':
    import argparse
    import json
    import os
    from subprocess import Popen
    import signal
    import time
    import shlex

    #process commandline Arguments
    parser = get_argument_parser()
    args = parser.parse_args()

    #dump argument to config.json
    generate_config(args)

    #run snakemake
    call_snakemake(args.RCAS_path, args.forcerun, args.cores)
