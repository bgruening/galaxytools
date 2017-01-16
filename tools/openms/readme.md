Galaxy wrapper for OpenMS
=========================

OpenMS is an open-source software C++ library for LC/MS data management and analyses.
It offers an infrastructure for the rapid development of mass spectrometry related software.
OpenMS is free software available under the three clause BSD license and runs under Windows, MacOSX and Linux.

More informations are available at:

 * https://github.com/OpenMS/OpenMS
 * http://open-ms.sourceforge.net


Installation
============

Galaxy should be able to automatically install the dependencies, i.e.
'package_openms_2_0' or 'package_qt_4_8' repository.

The wrappers are included in https://testtoolshed.g2.bx.psu.edu/view/bgruening/openms.


Generating OpenMS wrappers
==========================

 * install OpenMS (you can do this automatically through the Tool Shed)
 * create a folder called CTD
 * inside of your new installed openms/bin folder, execute the following command:
    
    ```bash
    for binary in `ls`; do ./$binary -write_ctd /PATH/TO/YOUR/CTD; done;
    ```
    
 * `MetaProSIP.ctd` includes a not supported character: To use it, search for `²` and replace it (e.g. with `^2`).

 * clone or install CTDopts

    ```bash
    git clone https://github.com/genericworkflownodes/CTDopts
    ```

 * add CTDopts to your `$PYTHONPATH`

    ```bash
    export PYTHONPATH=/home/user/CTDopts/
    ```

 * clone or install CTD2Galaxy

    ```bash
    git clone https://github.com/blankclemens/CTD2Galaxy.git
    ```
    
 * If you have CTDopts and CTD2Galaxy installed you are ready to generate Galaxy Tools from CTD definitions. Change the following command according to your needs, especially the `/PATH/TO` parts. 

    ```bash
    python generator.py \ 
    -i /PATH/TO/YOUR/CTD/*.ctd \
    -o ./PATH/TO/YOUR/WRAPPERS/ -t tool.conf \
    -d OpenMS -g openms \
    -b version log debug test no_progress threads \
     in_type executable myrimatch_executable \
     fido_executable fidocp_executable \
     omssa_executable pepnovo_executable \
     xtandem_executable param_model_directory \
     java_executable java_memory java_permgen \
    -f /PATH/TO/filetypes.txt -m /PATH/TO/macros.xml \
    -s PATH/TO/SKIP_TOOLS_FILES.txt
    ```


 * As last step you need to change manually the binary names of all external binaries you want to use in OpenMS. For example:

    ```
    sed -i '13 a\-fido_executable fido' wrappers/FidoAdapter.xml
    sed -i '13 a\-fidocp_executable fido_choose_parameters' wrappers/FidoAdapter.xml
    sed -i '13 a\-executable msgfplus.jar' wrappers/MSGFPlusAdapter.xml
    sed -i '13 a\-myrimatch_executable myrimatch' wrappers/MyriMatchAdapter.xml
    sed -i '13 a\-omssa_executable omssa' wrappers/OMSSAAdapter.xml
    sed -i '13 a\-pepnovo_executable pepnovo' wrappers/PepNovoAdapter.xml
    sed -i '13 a\-xtandem_executable xtandem' wrappers/XTandemAdapter.xml
    ```
    
 * For some tools, additional work has to be done. In `MSGFPlusAdapter.xml` the following is needed in the command section at the beginning (check your file to know what to copy where):
 
   ```
    <command><![CDATA[

    ## check input file type
    #set $in_type = $param_in.ext

    ## create the symlinks to set the proper file extension, since msgf uses them to choose how to handle the input files
    ln -s '$param_in' 'param_in.${in_type}' &&
    ln -s '$param_database' param_database.fasta &&
    ## find location of the MSGFPlus.jar file of the msgf_plus conda package
    MSGF_JAR=\$(msgf_plus -get_jar_path) &&

    MSGFPlusAdapter
    -executable \$MSGF_JAR
    #if $param_in:
      -in 'param_in.${in_type}'
    #end if
    #if $param_out:
      -out $param_out
    #end if
    #if $param_mzid_out:
      -mzid_out $param_mzid_out
    #end if
    #if $param_database:
      -database param_database.fasta
    #end if
    
    [...]
    ]]>
    ```

 * These tools have multiple outputs (number of inputs = number of outputs) which is not yet supported in Galaxy-stable and are therefore in `SKIP_TOOLS_FILES.txt`:
    * SeedListGenerator
    * SpecLibSearcher
    * MapAlignerIdentification
    * MapAlignerPoseClustering
    * MapAlignerSpectrum
    * MapAlignerRTTransformer
    
 * Additionally cause of lacking dependencies, the following adapters have been removed in `SKIP_TOOLS_FILES.txt` as well:
    * OMSSAAdapter
    * MyrimatchAdapter


Licence (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

