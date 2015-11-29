#!/usr/bin/env bash

for filepath in `git diff --name-only $TRAVIS_COMMIT_RANGE`; do
    inspect_dir=`dirname $filepath`
    shed_found=$(find $inspect_dir -maxdepth 1 -name '.shed.yml' | wc -l)

    if [ $shed_found -ne 0 ]; then
        # Here we can add very extensive testing
        # like starting a Docker container and testing against it
        planemo shed_lint --report_level warn --tools --fail_level error $shed_found
    fi
done
