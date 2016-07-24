import os
import shutil
import sys
import tempfile

FIJI_JAR_DIR = os.environ.get('FIJI_JAR_DIR', None)
FIJI_OSX_JAVA3D_DIR = os.environ.get('FIJI_OSX_JAVA3D_DIR', None)
FIJI_PLUGIN_DIR = os.environ.get('FIJI_PLUGIN_DIR', None)
FIJI_ROOT_DIR = os.environ.get('FIJI_ROOT_DIR', None)

BUFF_SIZE = 1048576


def cleanup_before_exit(tmp_dir):
    """
    Remove temporary files and directories prior to tool exit.
    """
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)


def get_base_cmd_bunwarpj(jvm_memory):
    if FIJI_JAR_DIR is not None and FIJI_PLUGIN_DIR is not None:
        if jvm_memory in [None, 'None']:
            jvm_memory_str = ''
        else:
            jvm_memory_str = '-Xmx%s' % jvm_memory
        # bunwarpj_base_cmd = "java %s -cp %s/ij-1.49k.jar:%s/bUnwarpJ_-2.6.1.jar bunwarpj.bUnwarpJ_" % \
        #    (jvm_memory_str, FIJI_JAR_DIR, FIJI_PLUGIN_DIR)
        bunwarpj_base_cmd = "bunwarpj %s" % jvm_memory_str
        return bunwarpj_base_cmd
    return None


def get_base_command_imagej2(memory_size=None, macro=None, jython_script=None):
    imagej2_executable = get_imagej2_executable()
    if imagej2_executable is None:
        return None
    cmd = '%s --ij2 --headless --debug' % imagej2_executable
    if memory_size is not None:
        memory_size_cmd = ' -DXms=%s -DXmx=%s' % (memory_size, memory_size)
        cmd += memory_size_cmd
    if macro is not None:
        cmd += ' --macro %s' % os.path.abspath(macro)
    if jython_script is not None:
        cmd += ' --jython %s' % os.path.abspath(jython_script)
    return cmd


def get_file_extension(image_format):
    """
    Return a valid bioformats file extension based on the received
    value of image_format(e.g., "gif" is returned as ".gif".
    """
    return '.%s' % image_format


def get_file_name_without_extension(file_path):
    """
    Eliminate the .ext from the received file name, assuming that
    the file name consists of only a single '.'.
    """
    if os.path.exists(file_path):
        path, name = os.path.split(file_path)
        name_items = name.split('.')
        return name_items[0]
    return None


def get_imagej2_executable():
    """
    Fiji names the ImageJ executable different names for different
    architectures, but our bioconda recipe allows us to do this.
    """
    return 'ImageJ'


def get_input_image_path(tmp_dir, input_file, image_format):
    """
    Bioformats uses file extensions (e.g., .job, .gif, etc)
    when reading and writing image files, so the Galaxy dataset
    naming convention of setting all file extensions as .dat
    must be handled.
    """
    image_path = get_temporary_image_path(tmp_dir, image_format)
    # Remove the file so we can create a symlink.
    os.remove(image_path)
    os.symlink(input_file, image_path)
    return image_path


def get_platform_info_dict():
    '''Return a dict with information about the current platform.'''
    platform_dict = {}
    sysname, nodename, release, version, machine = os.uname()
    platform_dict['os'] = sysname.lower()
    platform_dict['architecture'] = machine.lower()
    return platform_dict


def get_stderr_exception(tmp_err, tmp_stderr, tmp_out, tmp_stdout, include_stdout=False):
    tmp_stderr.close()
    """
    Return a stderr string of reasonable size.
    """
    # Get stderr, allowing for case where it's very large.
    tmp_stderr = open(tmp_err, 'rb')
    stderr_str = ''
    buffsize = BUFF_SIZE
    try:
        while True:
            stderr_str += tmp_stderr.read(buffsize)
            if not stderr_str or len(stderr_str) % buffsize != 0:
                break
    except OverflowError:
        pass
    tmp_stderr.close()
    if include_stdout:
        tmp_stdout = open(tmp_out, 'rb')
        stdout_str = ''
        buffsize = BUFF_SIZE
        try:
            while True:
                stdout_str += tmp_stdout.read(buffsize)
                if not stdout_str or len(stdout_str) % buffsize != 0:
                    break
        except OverflowError:
            pass
    tmp_stdout.close()
    if include_stdout:
        return 'STDOUT\n%s\n\nSTDERR\n%s\n' % (stdout_str, stderr_str)
    return stderr_str


def get_temp_dir(prefix='tmp-imagej-', dir=None):
    """
    Return a temporary directory.
    """
    return tempfile.mkdtemp(prefix=prefix, dir=dir)


def get_tempfilename(dir=None, suffix=None):
    """
    Return a temporary file name.
    """
    fd, name = tempfile.mkstemp(suffix=suffix, dir=dir)
    os.close(fd)
    return name


def get_temporary_image_path(tmp_dir, image_format):
    """
    Return the path to a temporary file with a valid image format
    file extension that can be used with bioformats.
    """
    file_extension = get_file_extension(image_format)
    return get_tempfilename(tmp_dir, file_extension)


def handle_none_type(val, val_type='float'):
    if val is None:
        return ' None'
    else:
        if val_type == 'float':
            return ' %.3f' % val
        elif val_type == 'int':
            return ' %d' % val
    return ' %s' % val


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)
