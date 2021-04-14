#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import ConfigParser


def manage_loc_file(path, data):
    """
    updates or inserts config lines like the following one:
    chembl_2011 ChemBL 2011 (~1.0M compounds)   schema  database_name   user    localhost       user_password
    """
    new_config_entry = "%s\t%s\t%s\tmoleculedb\t%s\tlocalhost\t%s" % (
        data["db_key"],
        data["description"],
        data["schema"],
        data["user"],
        data["password"],
    )
    lines = list()
    replace = False
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                lines.append(line.strip())
            else:
                db_key = line.strip().split("\t", 0)
                if db_key == data["db_key"].strip():
                    lines.append(new_config_entry)
                    replace = True
                else:
                    lines.append(line.strip())
    if not replace:
        lines.append(new_config_entry)
    with open(path, "w") as handle:
        handle.write("\n".join(lines))


loc_file = os.path.join(
    sys.argv[5], "tool-data", "cheminformatics_chemical_searches_db.loc"
)

if sys.argv[1] == "--non-admin-user":
    error_text = (
        "The user %s have no permissions to insert data in the database, please contact the administrator or check your settings in universe_wsgi.ini."
        % sys.argv[2]
    )
    open(sys.argv[3], "w").write(error_text)
    sys.stdout.write(error_text)
    sys.exit()

config = ConfigParser.RawConfigParser()
config.read(os.path.join(sys.argv[5], "universe_wsgi.ini"))
database_connection = config.get("app:main", "database_connection")
user, password = database_connection.split("//")[-1].split("@")[0].split(":")

outfile = open(sys.argv[4], "w")
processes = sys.argv[3]

exec_dir = os.path.split(sys.argv[0])[0]

tempdir = os.path.abspath(tempfile.mkdtemp())
executeable = os.path.join(tempdir, "GenericMoleculeParser.py")
db_definition = os.path.join(tempdir, "GenericMoleculeDB.py")
utils_path = os.path.join(tempdir, "cheminfolib.py")
shutil.copy(os.path.join(exec_dir, "__init__.py"), os.path.join(tempdir, "__init__.py"))
shutil.copy(os.path.join(exec_dir, "GenericMoleculeParser.py"), executeable)
shutil.copy(os.path.join(exec_dir, "GenericMoleculeDB.py"), db_definition)
shutil.copy(os.path.join(exec_dir, "cheminfolib.py"), utils_path)


"""
We need to write a python config file, that can be imported into the parser.
It's not easy to give a declarative base class a scheme on runtime only on compile-time.
"""
gfile = open(os.path.join(tempdir, "generic.py"), "w")
gfile.write('host = "localhost"\n')
gfile.write('username = "%s"\n' % user)
gfile.write('password = "%s"\n' % password)
gfile.write('database = "moleculedb"\n')
gfile.write('schema = "%s"\n' % sys.argv[2])
gfile.close()

# Update the loc config file for all databases
config_data = {
    "user": user,
    "password": password,
    "db_key": sys.argv[2],
    "schema": sys.argv[2],
    "description": sys.argv[6],
}
manage_loc_file(loc_file, config_data)

try:
    abs_exec_path = os.path.abspath(executeable)
    args = [
        abs_exec_path,
        "--processes",
        processes,
        "--molecule-file",
        sys.argv[1],
        "--filetype",
        sys.argv[7],
    ]
    child = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except Exception as err:
    sys.stderr.write(
        "Return error code %i from command:\n%s\n" % (return_code, " ".join(args))
    )
    sys.exit(1)

stdout, stderr = child.communicate()
return_code = child.returncode

if return_code:
    sys.stdout.write(stdout)
    sys.stderr.write(stderr)
    sys.stderr.write(
        "Return error code %i from command:\n%s\n" % (return_code, " ".join(args))
    )
else:
    sys.stdout.write(stdout)
    outfile.write(
        "The file %s was sucessfully inserted into the database.\nPlease remember to check your .loc file under ./tool-data/ and restart galaxy.\n"
        % (sys.argv[2])
    )
    if stderr:
        outfile.write(
            "The following structures failed to import to the database:\n%s" % (stderr)
        )
    outfile.close()

sys.exit()
