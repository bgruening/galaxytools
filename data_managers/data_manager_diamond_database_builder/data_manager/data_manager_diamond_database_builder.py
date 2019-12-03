#!/usr/bin/env python
import json
import sys
import os
import tempfile
import shutil
import optparse
import urllib2
import subprocess
from ftplib import FTP
import tarfile
import zipfile
import gzip
import bz2

CHUNK_SIZE = 2**20  # 1mb


def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def _get_files_in_ftp_path(ftp, path):
    path_contents = []
    ftp.retrlines('MLSD %s' % (path), path_contents.append)
    return [line.split(';')[-1].lstrip() for line in path_contents]


def _get_stream_readers_for_tar(file_obj, tmp_dir):
    fasta_tar = tarfile.open(fileobj=file_obj, mode='r:*')
    return [fasta_tar.extractfile(member) for member in fasta_tar.getmembers()]


def _get_stream_readers_for_zip(file_obj, tmp_dir):
    fasta_zip = zipfile.ZipFile(file_obj, 'r')
    rval = []
    for member in fasta_zip.namelist():
        fasta_zip.extract(member, tmp_dir)
        rval.append(open(os.path.join(tmp_dir, member), 'rb'))
    return rval


def _get_stream_readers_for_gzip(file_obj, tmp_dir):
    return [gzip.GzipFile(fileobj=file_obj, mode='rb')]


def _get_stream_readers_for_bz2(file_obj, tmp_dir):
    return [bz2.BZ2File(file_obj.name, 'rb')]


def download_from_ncbi(data_manager_dict, params, target_directory, database_id, database_name):
    NCBI_FTP_SERVER = 'ftp.ncbi.nlm.nih.gov'
    NCBI_DOWNLOAD_PATH = '/blast/db/FASTA/'
    COMPRESSED_EXTENSIONS = [('.tar.gz', _get_stream_readers_for_tar), ('.tar.bz2', _get_stream_readers_for_tar), ('.zip', _get_stream_readers_for_zip), ('.gz', _get_stream_readers_for_gzip), ('.bz2', _get_stream_readers_for_bz2)]

    ncbi_identifier = params['param_dict']['reference_source']['requested_identifier']
    ftp = FTP(NCBI_FTP_SERVER)
    ftp.login()

    path_contents = _get_files_in_ftp_path(ftp, NCBI_DOWNLOAD_PATH)

    ncbi_file_name = None
    get_stream_reader = None
    ext = None
    for ext, get_stream_reader in COMPRESSED_EXTENSIONS:
        if "%s%s" % (ncbi_identifier, ext) in path_contents:
            ncbi_file_name = "%s%s%s" % (NCBI_DOWNLOAD_PATH, ncbi_identifier, ext)
            break

    if not ncbi_file_name:
        raise Exception('Unable to determine filename for NCBI database for %s: %s' % (ncbi_identifier, path_contents))

    tmp_dir = tempfile.mkdtemp(prefix='tmp-data-manager-ncbi-')
    ncbi_fasta_filename = os.path.join(tmp_dir, "%s%s" % (ncbi_identifier, ext))

    fasta_base_filename = "%s.fa" % database_id
    fasta_filename = os.path.join(target_directory, fasta_base_filename)
    fasta_writer = open(fasta_filename, 'wb+')

    tmp_extract_dir = os.path.join(tmp_dir, 'extracted_fasta')
    os.mkdir(tmp_extract_dir)

    tmp_fasta = open(ncbi_fasta_filename, 'wb+')

    ftp.retrbinary('RETR %s' % ncbi_file_name, tmp_fasta.write)

    tmp_fasta.flush()
    tmp_fasta.seek(0)

    fasta_readers = get_stream_reader(tmp_fasta, tmp_extract_dir)

    data_table_entry = _stream_fasta_to_file(fasta_readers, target_directory, database_id, database_name, params)
    _add_data_table_entry(data_manager_dict, data_table_entry)

    for fasta_reader in fasta_readers:
        fasta_reader.close()
    tmp_fasta.close()
    cleanup_before_exit(tmp_dir)


def download_from_url(data_manager_dict, params, target_directory, database_id, database_name):
    # TODO: we should automatically do decompression here
    urls = filter(bool, map(lambda x: x.strip(), params['param_dict']['reference_source']['user_url'].split('\n')))
    fasta_reader = [urllib2.urlopen(url) for url in urls]

    data_table_entry = _stream_fasta_to_file(fasta_reader, target_directory, database_id, database_name, params)
    _add_data_table_entry(data_manager_dict, data_table_entry)


def download_from_history(data_manager_dict, params, target_directory, database_id, database_name):
    # TODO: allow multiple FASTA input files
    input_filename = params['param_dict']['reference_source']['input_fasta']
    if isinstance(input_filename, list):
        fasta_reader = [open(filename, 'rb') for filename in input_filename]
    else:
        fasta_reader = open(input_filename)

    data_table_entry = _stream_fasta_to_file(fasta_reader, target_directory, database_id, database_name, params)
    _add_data_table_entry(data_manager_dict, data_table_entry)


def copy_from_directory(data_manager_dict, params, target_directory, database_id, database_name):
    input_filename = params['param_dict']['reference_source']['fasta_filename']
    create_symlink = params['param_dict']['reference_source']['create_symlink'] == 'create_symlink'
    if create_symlink:
        data_table_entry = _create_symlink(input_filename, target_directory, database_id, database_name)
    else:
        if isinstance(input_filename, list):
            fasta_reader = [open(filename, 'rb') for filename in input_filename]
        else:
            fasta_reader = open(input_filename)
        data_table_entry = _stream_fasta_to_file(fasta_reader, target_directory, database_id, database_name, params)
    _add_data_table_entry(data_manager_dict, data_table_entry)


def _add_data_table_entry(data_manager_dict, data_table_entry):
    data_manager_dict['data_tables'] = data_manager_dict.get('data_tables', {})
    data_manager_dict['data_tables']['diamond_database'] = data_manager_dict['data_tables'].get('diamond_database', [])
    data_manager_dict['data_tables']['diamond_database'].append(data_table_entry)
    return data_manager_dict


def _stream_fasta_to_file(fasta_stream, target_directory, database_id, database_name, params, close_stream=True):
    fasta_base_filename = "%s.fa" % database_id
    fasta_filename = os.path.join(target_directory, fasta_base_filename)

    temp_fasta = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
    temp_fasta.close()
    fasta_writer = open(temp_fasta.name, 'wb+')

    if isinstance(fasta_stream, list) and len(fasta_stream) == 1:
        fasta_stream = fasta_stream[0]

    if isinstance(fasta_stream, list):
        last_char = None
        for fh in fasta_stream:
            if last_char not in [None, '\n', '\r']:
                fasta_writer.write('\n')
            while True:
                data = fh.read(CHUNK_SIZE)
                if data:
                    fasta_writer.write(data)
                    last_char = data[-1]
                else:
                    break
            if close_stream:
                fh.close()
    else:
        while True:
            data = fasta_stream.read(CHUNK_SIZE)
            if data:
                fasta_writer.write(data)
            else:
                break
        if close_stream:
            fasta_stream.close()

    fasta_writer.close()

    args = ['diamond', 'makedb', '--in', temp_fasta.name, '--db', fasta_filename]

    tmp_stderr = tempfile.NamedTemporaryFile(prefix="tmp-data-manager-diamond-database-builder-stderr")
    proc = subprocess.Popen(args=args, shell=False, cwd=target_directory, stderr=tmp_stderr.fileno())
    return_code = proc.wait()
    if return_code:
        tmp_stderr.flush()
        tmp_stderr.seek(0)
        print >> sys.stderr, "Error building diamond database:"
        while True:
            chunk = tmp_stderr.read(CHUNK_SIZE)
            if not chunk:
                break
            sys.stderr.write(chunk)
        sys.exit(return_code)
    tmp_stderr.close()
    os.remove(temp_fasta.name)
    return dict(value=database_id, name=database_name, db_path="%s.dmnd" % fasta_base_filename)


def _create_symlink(input_filename, target_directory, database_id, database_name):
    fasta_base_filename = "%s.fa" % database_id
    fasta_filename = os.path.join(target_directory, fasta_base_filename)
    os.symlink(input_filename, fasta_filename)
    return dict(value=database_id, name=database_name, db_path=fasta_base_filename)


REFERENCE_SOURCE_TO_DOWNLOAD = dict(ncbi=download_from_ncbi, url=download_from_url, history=download_from_history, directory=copy_from_directory)


def main():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dbkey_description', dest='dbkey_description', action='store', type="string", default=None, help='dbkey_description')
    (options, args) = parser.parse_args()

    filename = args[0]

    params = json.loads(open(filename).read())
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory)
    data_manager_dict = {}

    database_id = params['param_dict']['database_id']
    database_name = params['param_dict']['database_name']

    # Fetch the FASTA
    REFERENCE_SOURCE_TO_DOWNLOAD[params['param_dict']['reference_source']['reference_source_selector']](data_manager_dict, params, target_directory, database_id, database_name)

    # save info to json file
    open(filename, 'wb').write(json.dumps(data_manager_dict))


if __name__ == "__main__":
    main()
