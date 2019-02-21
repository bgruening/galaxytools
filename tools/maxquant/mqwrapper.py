"""
Run MaxQuant on a modified mqpar.xml.
Use maxquant conda package.
TODO: - extend outputs
      (- add support for protein groups)
      - use file.displayname to match raw file names
      - create scripts to read modifications.xml and edit mqwrapper.xml

Authors: Damian Glaetzer <d.glaetzer@mailbox.org>

based on the maxquant galaxy tool by John Chilton:
https://github.com/galaxyproteomics/tools-galaxyp/tree/master/tools/maxquant
"""

import os
import argparse
import subprocess
import shutil
import mqparam

# build parser
parser = argparse.ArgumentParser()
arguments = ["--raw_files", "--fasta_files", "--fixed_mods",
             "--var_mods", "--proteases", "--exp_design",
             "--missed_cleavages", "--min_unique_pep", "--mqpar_in",
             "--num_threads", "--tool_dir", "--output_all", "--mqpar_out",
             "--raw_file_names", "--mzTab", "--light_mods", "--medium_mods",
             "--heavy_mods"]

flags = ("--calc_peak_properties", "--silac", "--write_mztab")


txt_output = ("evidence", "msms", "parameters",
              "peptides", "proteinGroups", "allPeptides",
              "libraryMatch", "matchedFeatures",
              "modificationSpecificPeptides", "ms3Scans",
              "msmsScans", "mzRange", "peptideSection", "summary")

for el in txt_output:
    arguments.append('--' + el)

for arg in arguments:
    parser.add_argument(arg)
for flag in flags:
    parser.add_argument(flag, action="store_true")

args = vars(parser.parse_args())

# link rawfile datasets to .raw names for maxquant to accept them
files = args['raw_files'].split(',')
filenames = args['raw_file_names'].split(',')
for f, l in zip(files, filenames):
    os.link(f, l)

# arguments for mqparam
simple_args = ('missed_cleavages', 'min_unique_pep',
               'num_threads', 'calc_peak_properties',
               'write_mztab')

list_args = ('fixed_mods', 'var_mods', 'proteases')

# build mqpar.xml
template_name = os.path.join(args['tool_dir'], 'template.xml')
mqpar_in = args['mqpar_in'] if args['mqpar_in'] else template_name
mqpar_temp = os.path.join(os.getcwd(), 'mqpar.xml')
mqpar_out = args['mqpar_out'] if args['mqpar_out'] else mqpar_temp

m = mqparam.MQParam(mqpar_out, mqpar_in, args['exp_design'])
m.add_rawfiles([os.path.join(os.getcwd(), name) for name in filenames])
m.add_fasta_files(args['fasta_files'].split(','))

for e in simple_args:
    if args[e]:
        m.set_simple_param(e, args[e])

for e in list_args:
    if args[e]:
        m.set_list_params(e, args[e].split(','))

m.set_silac(args['light_mods'].split(',') if args['light_mods'] else None,
            args['medium_mods'].split(',') if args['medium_mods'] else None,
            args['heavy_mods'].split(',') if args['heavy_mods'] else None)

m.write()

# build and run MaxQuant command
cmd = ['maxquant', args['mqpar_out']]

subprocess.run(cmd, check=True, cwd='./')

# copy results to galaxy database
for el in txt_output:
    destination = args[el]
    source = os.path.join(os.getcwd(), "combined", "txt", "{}.txt".format(el))
    if destination and os.path.isfile(source):
        shutil.copy(source, destination)

if args['mzTab']:
    source = os.path.join(os.getcwd(), "combined", "txt", "mzTab.mzTab")
    if os.path.isfile(source):
        shutil.copy(source, args['mzTab'])

if args['output_all']:
    subprocess.run(('tar', '-zcf', args['output_all'], './combined/txt/'))
