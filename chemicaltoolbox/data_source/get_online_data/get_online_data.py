import os
import urllib.request
import gzip, tempfile
import zipfile
import subprocess
import shutil
import argparse
from io import BytesIO

def unescape(cond_text):
    # Unescape if input has been escaped
    mapped_chars = { '>' :'__gt__', 
                 '<' :'__lt__', 
                 "'" :'__sq__',
                 '"' :'__dq__',
                 '[' :'__ob__',
                 ']' :'__cb__',
                 '{' :'__oc__',
                 '}' :'__cc__',
                 '@' : '__at__',
                 '\n' : '__cn__',
                 '\r' : '__cr__',
                 '\t' : '__tc__'
                 }
    for key, value in mapped_chars.items():
        cond_text = cond_text.replace( value, key )
    return cond_text

def get_files(options):
    urls = unescape(options.url)
    with open(options.out, 'wb+') as out:
        if options.whitelist:
            allowed_extensions = [ext.strip() for ext in unescape(options.whitelist).split('\n')]
        else:
            allowed_extensions = ['.sdf', '.smi', '.inchi', '.mol']

        for url in urls.split('\n'):
            request = urllib.request.Request(url)
            response = urllib.request.urlopen(request)
            resp_read = response.read()
            if resp_read[:2] == b'\x1f\x8b':  # test magic number for gzipped files
                response = urllib.request.urlopen(request)
                out.write(gzip.decompress(resp_read))
            elif resp_read[:2] == b'PK':  # test magic number for zipped files
                temp = tempfile.NamedTemporaryFile(delete=False)
                temp.close()
                zf = zipfile.ZipFile(BytesIO(resp_read), allowZip64=True)
                tmpdir = tempfile.mkdtemp()

                for filename in zf.namelist():
                    zf.extractall(tmpdir)

                os.remove(temp.name)
                molfiles = []
                for root, dirs, files in os.walk(tmpdir):
                    for filename in files:
                        if os.path.splitext(filename)[-1].lower() in allowed_extensions or allowed_extensions == []:
                            mfile = os.path.join(root, filename)
                            shutil.copyfileobj(open(mfile, 'rb'), out)
                shutil.rmtree( tmpdir )
                zf.close()
            else:
                out.write(resp_read)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Download compressed files and extract files of with chosen extensions
    """)
    parser.add_argument('--url', dest='url', help='URL')
    parser.add_argument('--whitelist', dest='whitelist', default=None, help='whitelist')
    parser.add_argument('--out', dest='out', help='output')
    
    options = parser.parse_args()
    get_files(options)