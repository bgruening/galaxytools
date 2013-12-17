#!/bin/sh

##
## Try to build GNU coreutils from source
##
URL=http://ftp.gnu.org/gnu/coreutils/coreutils-8.16.tar.xz
FILE=$(basename "$URL") || exit 1
DIR=$(basename "$FILE" .tar.xz) || exit 1

PREFIX=$PWD

mkdir -p "$PREFIX/usr" || exit 1

rm -rf "$FILE" "$DIR"

echo "Downloading coreutils from $URL..."
wget -q "$URL" || exit 1
echo "Extracting source files..."
tar -xJf "$FILE" || exit 1

cd "$DIR" || exit 1
./configure --prefix "$PREFIX/usr" || exit 1
make || exit 1
make install || exit 1
cd ".." || exit 1

echo "

Done!
compiled executables are in $PREFIX/usr/bin
Either add them to the Galaxy PATH or copy them to a PATH directory

"
