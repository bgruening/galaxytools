#!/bin/sh

##
## Try to build GNU SED from source (+ patch it)
##
URL=ftp://ftp.gnu.org/gnu/sed/sed-4.2.tar.gz
FILE=$(basename "$URL") || exit 1
DIR=$(basename "$FILE" .tar.gz) || exit 1

PREFIX=$PWD

mkdir -p "$PREFIX/usr" || exit 1

rm -rf "$FILE" "$DIR"

echo "Downloading SED from $URL..."
wget -q "$URL" || exit 1
echo "Extracting source files..."
tar -xzf "$FILE" || exit 1

cd "$DIR" || exit 1
patch -p1 < "$PREFIX/patches/sed-4.2-sandbox.patch" || exit 1
./configure --prefix "$PREFIX/usr" || exit 1
make || exit 1
make install || exit 1
cd ".." || exit 1

echo "

Done!
compiled executables are in $PREFIX/usr/bin
Either add them to the Galaxy PATH or copy them to a PATH directory

"

