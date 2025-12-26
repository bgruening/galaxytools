To run Galaxy tests, the following environment variables must be set:

BIOBLEND_GALAXY_API_KEY = <web API key of a running Galaxy instance>
BIOBLEND_GALAXY_URL = <root url for a running Galaxy instance like http://localhost:8080>

To run ToolShed tests, the following environment variable must be set:

BIOBLEND_TOOLSHED_URL = <root url for a running ToolShed instance like https://testtoolshed.g2.bx.psu.edu/>

If you wish to run the entire suite, set all of the above.  The integration
tests can subsequently be run by invoking `pytest` from the command line.
pytest should be invoked from the project root folder, and not the tests
child folder, since the test data is resolved relative to the bioblend folder.
