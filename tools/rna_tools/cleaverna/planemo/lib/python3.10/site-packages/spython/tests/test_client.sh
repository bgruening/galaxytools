#!/bin/bash

# Include help functions
. helpers.sh

echo
echo "************** START: test_client.sh **********************"

# Create temporary testing directory
echo "Creating temporary directory to work in."
tmpdir=$(mktemp -d)
output=$(mktemp ${tmpdir:-/tmp}/spython_test.XXXXXX)
here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "Testing help commands..."

# Test help for all commands
for command in recipe shell;
    do
    runTest 0 $output spython $command --help
done

echo "#### Testing recipe auto generation"
runTest 1 $output spython recipe $here/testdata/Dockerfile | grep "FROM"
runTest 0 $output spython recipe $here/testdata/Dockerfile | grep "%post"
runTest 1 $output spython recipe $here/testdata/Singularity | grep "%post"
runTest 0 $output spython recipe $here/testdata/Singularity | grep "FROM"

echo "#### Testing recipe targeted generation"
runTest 0 $output spython recipe --writer docker $here/testdata/Dockerfile | grep "FROM"
runTest 1 $output spython recipe --writer docker $here/testdata/Dockerfile | grep "%post"
runTest 0 $output spython recipe --writer singularity $here/testdata/Singularity | grep "%post"
runTest 1 $output spython recipe --writer singularity $here/testdata/Singularity | grep "FROM"

echo "#### Testing recipe file generation"
outfile=$(mktemp ${tmpdir:-/tmp}/spython_recipe.XXXXXX)
runTest 0 $output spython recipe $here/testdata/Dockerfile $outfile
runTest 0 $output test -f "$outfile"
runTest 0 $output cat $outfile | grep "%post"
rm $outfile

echo "#### Testing recipe json export"
runTest 0 $output spython recipe --json $here/testdata/Dockerfile | grep "ports"
runTest 0 $output spython recipe $here/testdata/Dockerfile $outfile
runTest 0 $output test -f "$outfile"
runTest 0 $output cat $outfile | grep "%post"

# Force is false, should fail
echo "#### Testing recipe json export, writing to file"
runTest 0 $output spython recipe --json $here/testdata/Dockerfile $outfile
runTest 0 $output spython recipe --force --json $here/testdata/Dockerfile $outfile
runTest 0 $output test -f "$outfile"
runTest 0 $output cat $outfile | grep "ports"

echo "Finish testing basic client"
rm -rf ${tmpdir}
