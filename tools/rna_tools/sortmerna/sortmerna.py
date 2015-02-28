#!/usr/bin/env python

"""
Runs SortMeRNA
"""

import subprocess
import optparse
import shlex


def main():
    """Parse the command line, exectutes SortMeRNA and buildtrie if neeeded."""
    #TODO: Put all SortMeRNA options in the command-line parser
    parser = optparse.OptionParser()
    parser.add_option('--sortmerna', dest='sortmerna_cmd', help='')
    parser.add_option('--buildtrie', dest='buildtrie',
                      default=False, action='store_true', help='')
    (options, args) = parser.parse_args()
    if not args:
        raise Exception('Please provide at least one database')

    if options.buildtrie:
        buildtrie = 'buildtrie'
        for database in args:
            run_buildtrie([buildtrie, '--db', database])

    if options.sortmerna_cmd:
        sortmerna = 'sortmerna'
        run_sortmerna([sortmerna] +
                      shlex.split(options.sortmerna_cmd) +
                      ['-m', '262144', '-n', str(len(args)), '--db'] +
                      args)


def run_buildtrie(cmd):
    """Run the BuildTrie program."""
    try:
        stdout_arg = subprocess.PIPE
        stderr_arg = subprocess.PIPE
        child_process = subprocess.Popen(args=" ".join(cmd), shell=True,
                                         stdin=None, stdout=stdout_arg,
                                         stderr=stderr_arg)
        stdout_str, stderr_str = child_process.communicate()
        return_code = child_process.returncode
        if return_code is not 0:
            raise Exception(stderr_str)

    except Exception, error:
        raise Exception('Error while running Buildtrie:\n' +
                        '\n'.join([str(error), stdout_str, stderr_str]))


def run_sortmerna(cmd):
    """Run the SortMeRNA program."""
    try:
        stdout_arg = subprocess.PIPE
        stderr_arg = subprocess.PIPE
        child_process = subprocess.Popen(args=" ".join(cmd), shell=True,
                                         stdin=None, stdout=stdout_arg,
                                         stderr=stderr_arg)
        stdout_str, stderr_str = child_process.communicate()
        return_code = child_process.returncode
        if return_code is not 0:
            raise Exception(stderr_str)
    except Exception, error:
        raise Exception('Error while running SortMeRNA:\n' +
                        '\n'.join([str(error), stdout_str, stderr_str]))


if __name__ == "__main__":
    main()
