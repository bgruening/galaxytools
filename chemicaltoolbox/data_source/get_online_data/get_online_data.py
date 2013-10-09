#!/usr/bin/env python

__author__ = 'Bjoern Gruening'
__version__ = '0.1'
__date__ = '2012'
__license__ = 'GLP3+'

import os, sys
import urllib2
import gzip, tempfile
import zipfile
import subprocess
import shutil

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

urls = unescape(sys.argv[1])
out = open(sys.argv[2], 'wb')

if len(sys.argv) > 3:
    allowed_extensions = [ ext.strip() for ext in unescape(sys.argv[3]).split('\n') ]
else:
    allowed_extensions = ['.sdf', '.smi', '.inchi', '.mol']

for url in urls.split('\n'):
    url = url.strip()
    request = urllib2.Request( url )
    request.add_header('Accept-encoding', 'gzip')
    request.add_header('Accept-encoding', 'gz')
    response = urllib2.urlopen( request )

    if response.info().get('Content-Encoding') in ['gz','gzip'] or os.path.splitext(url)[-1] in ['.gz','.gzip']:
        temp = tempfile.NamedTemporaryFile( delete=False )
        temp.write( response.read() )
        temp.close()
        zipfile = gzip.open(temp.name, 'rb')
        out.write( zipfile.read() )
        os.remove(temp.name)
    elif response.info().get('Content-Encoding') in ['zip'] or os.path.splitext(url)[-1] in ['.zip']:
        temp = tempfile.NamedTemporaryFile(delete=False)
        temp.close()
        with open(temp.name, 'wb') as fp:
            shutil.copyfileobj(response, fp)

        zf = zipfile.ZipFile(temp.name, allowZip64=True)
        tmpdir = tempfile.mkdtemp( )

        for filename in zf.namelist():
            zf.extractall( tmpdir )

        os.remove( temp.name )
        molfiles = []
        for root, dirs, files in os.walk(tmpdir):
            for filename in files:
                if os.path.splitext(filename)[-1].lower() in allowed_extensions or allowed_extensions == []:
                    mfile = os.path.join( root, filename)
                    molfiles.append( mfile )

        for filename in molfiles:
            shutil.copyfileobj(open(filename, 'rb'), out)
        shutil.rmtree( tmpdir )
        zf.close()
    elif response.info().get('Content-Encoding') == 'rar' or os.path.splitext(url)[-1] in ['.rar']:
        temp = tempfile.NamedTemporaryFile(delete=False)
        temp.close()
        with open(temp.name, 'wb') as fp:
            shutil.copyfileobj(response, fp)
        cmd = subprocess.Popen('unrar p -inul %s' % temp.name, stdout=out, shell=True)
        os.remove( temp.name )
    else:
        out.write( response.read() )
out.close()
