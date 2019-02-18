"""
Create a project-specific MaxQuant parameter file.

TODO: make FASTA-template dynamic

Authors: Damian Glaetzer <d.glaetzer@mailbox.org>
         Franziska Elsaesser <fels@leute.server.de>
"""

import os
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom
from itertools import zip_longest


class MQParam:
    """Represents a mqpar.xml and provides methods to modify
    some of its parameters. Useful to either create a working
    mqpar.xml from a template or to just modify the file paths
    in an already complete mqpar-file.
    """

    fasta_template = """<FastaFileInfo>
    <fastaFilePath></fastaFilePath>
    <identifierParseRule></identifierParseRule>
    <descriptionParseRule></descriptionParseRule>
    <taxonomyParseRule></taxonomyParseRule>
    <variationParseRule></variationParseRule>
    <modificationParseRule></modificationParseRule>
    <taxonomyId></taxonomyId>
    </FastaFileInfo>"""

    version = '1.6.3.4'
    # map simple params to their node in the xml tree
    simple_params = {'missed_cleavages' :
                     '.parameterGroups/parameterGroup/maxMissedCleavages',
                     'min_unique_pep' : '.minUniquePeptides',
                     'num_threads' : 'numThreads',
                     'calc_peak_properties' : '.calcPeakProperties',
                     'write_mztab' : 'writeMzTab'}

    def __init__(self, mqpar_out, mqpar_in, exp_design):
        """Initialize MQParam class. mqpar_in can either be a template
        or a already suitable mqpar file.
        >>> t = MQParam("test", './template.xml', None)
        >>> t.root.tag
        'MaxQuantParams'
        >>> (t.root.find('maxQuantVersion')).text
        '1.6.3.4'
        """

        self.orig_mqpar = mqpar_in
        self.exp_design = exp_design
        self.mqpar_out = mqpar_out
        self.root = ET.parse(mqpar_in).getroot()
        v = self.root.find('maxQuantVersion').text
        if v != self.version:
            raise Exception("This tool is designed for MaxQuant "
                            + "{}.\n".format(self.version)
                            + "Your mqpar.xml is from MaxQuant {}".format(v))

    @staticmethod
    def _add_child(el, name, text, attrib=None):
        """Add a child element to an element.

        >>> t = MQParam("test", './template.xml', None)
        >>> MQParam._add_child(t.root, "test", "test")
        >>> t.root.find('test').text == "test"
        True
        """

        child = ET.SubElement(el, name, attrib=attrib if attrib else {})
        child.text = text

    @staticmethod
    def get_orig_filename(filename):
        """Search for the original file name in a raw file.
        Not the best solution for mapping raw files uploaded in
        galaxy to a locally created experimental design template."""

        # 1st group: file name without extension
        pattern = re.compile(r'[\\/]([^\./\\]+)\.raw')
        with open(filename, mode='rb') as f:
            # only look at first four lines
            for i in range(0, 4):
                line = f.readline()
                temp = []
                # ignore 0x00 bytes
                for c in line:
                    if c:
                        temp.append(c)
                ascii_line = bytes(temp).decode('ascii', 'ignore')
                match = re.search(pattern, ascii_line)
                if match:
                    return match.group(1)

        # raise error if file name could not be found
        raise Exception('Unable to find file name in {}'.format(filename))

    def add_rawfiles(self, rawfiles):
        """Add a list of raw files to the mqpar.xml.
        If experimental design template was specified,
        modify other parameters accordingly.
        The raw file must be specified as absolute paths
        for maxquant to find them.
        >>> t1 = MQParam("test", './template.xml', None)
        >>> t1.add_rawfiles(('test1', ))
        >>> t1.root.find("filePaths")[0].text
        'test1'
        >>> t1.root.find("fractions")[0].text
        '32767'
        >>> len(t1.root.find("fractions"))
        1
        >>> t2 = MQParam("test", './template.xml', \
                         './test-data/exp_design_test.txt')
        >>> t2.add_rawfiles(('test-data/QEplus021874.raw', \
                             'test-data/QEplus021876.raw'))
        >>> len(t2.root.find("filePaths"))
        2
        >>> t2.root.find("filePaths")[1].text
        'test-data/QEplus021876.raw'
        >>> t2.root.find("experiments")[1].text
        '2'
        >>> t2.root.find("fractions")[0].text
        '3'
        """

        # choose experimental design template
        if not self.exp_design:  # no experimentalDesignTemplate.txt given
            names = rawfiles
            fracs = ('32767',) * len(rawfiles)
            exps = [os.path.split(r)[1] for r in rawfiles]
            PTMs = ['False'] * len(rawfiles)
        else:  # experimentalDesignTemplate.txt given
            with open(self.exp_design) as design_file:
                design = {}
                index = []
                for line in design_file:
                    if design:
                        for e, i in zip_longest(line.strip().split('\t'), index):
                            design[i].append(e)
                    else:
                        for i in line.strip().split('\t'):
                            design[i] = []
                            index.append(i)

            # map rawfiles to names in exp. design template
            names = []
            names_to_paths = {}
            for r in rawfiles:
                names_to_paths[os.path.splitext(os.path.basename(r))[0]] = r
            for name in design['Name']:
                names.append(names_to_paths[name] if name in names_to_paths
                             else None)

            PTMslist = design['PTM']
            PTMs = []
            for item in PTMslist:
                if item == 'True':
                    PTMs.append(item)
                else:
                    PTMs.append('False')

            exps = design['Experiment']
            fracs = design['Fraction']

        # These parent nodes will get a child appended for each RAW file
        nodenames = ('filePaths', 'experiments', 'fractions',
                     'ptms', 'paramGroupIndices', 'referenceChannel')

        # Get parent nodes from document
        nodes = dict()
        for nodename in nodenames:
            node = self.root.find(nodename)
            if node is None:
                raise ValueError('Element {} not found in XML document'
                                 .format(nodename))
            nodes[nodename] = node
            node.clear()
            node.tag = nodename

        # Append sub-elements to nodes (one per file)
        for i in range(0, len(names)):
            if names[i]:
                MQParam._add_child(nodes['filePaths'], 'string', names[i])
                MQParam._add_child(nodes['experiments'], 'string', exps[i])
                MQParam._add_child(nodes['fractions'], 'short', fracs[i])
                MQParam._add_child(nodes['ptms'], 'boolean', PTMs[i])
                MQParam._add_child(nodes['paramGroupIndices'], 'int', '0')
                MQParam._add_child(nodes['referenceChannel'], 'string', '')

    def add_fasta_files(self, files):
        """Add fasta file groups.
        >>> t = MQParam('test', './template.xml', None)
        >>> t.add_fasta_files(('test1', 'test2'))
        >>> len(t.root.find('fastaFiles'))
        2
        >>> t.root.find('fastaFiles')[0].find("fastaFilePath").text
        'test1'
        """
        fasta_node = self.root.find("fastaFiles")
        fasta_node.clear()
        fasta_node.tag = "fastaFiles"

        identifier = r'>([^\s]*)'
        description = r'>(.*)'
        for index in range(len(files)):
            filepath = '<fastaFilePath>' + files[index]
            fasta = self.fasta_template.replace('<fastaFilePath>', filepath)
            fasta = fasta.replace('<identifierParseRule>',
                                  '<identifierParseRule>' + identifier)
            fasta = fasta.replace('<descriptionParseRule>',
                                  '<descriptionParseRule>' + description)
            ff_node = self.root.find('.fastaFiles')
            fastaentry = ET.fromstring(fasta)
            ff_node.append(fastaentry)

    def set_simple_param(self, key, value):
        """Set a simple parameter.
        >>> t = MQParam(None, './template.xml', None)
        >>> t.set_simple_param('min_unique_pep', 4)
        >>> t.root.find('.minUniquePeptides').text
        '4'
        """

        if key in self.simple_params:
            node = self.root.find(self.simple_params[key])
            node.text = str(value)
        else:
            raise ValueError("Parameter not found.")

    def set_silac(self, light_mods, medium_mods, heavy_mods):
        """Set label modifications.
        >>> t1 = MQParam('test', './template.xml', None)
        >>> t1.set_silac(None, ('test1', 'test2'), None)
        >>> t1.root.find('.parameterGroups/parameterGroup/maxLabeledAa').text
        '2'
        >>> t1.root.find('.parameterGroups/parameterGroup/multiplicity').text
        '3'
        >>> t1.root.find('.parameterGroups/parameterGroup/labelMods')[1].text
        'test1;test2'
        >>> t1.root.find('.parameterGroups/parameterGroup/labelMods')[2].text
        ''
        """
        multiplicity = 3 if medium_mods else 2 if heavy_mods else 1
        max_label = str(max(len(light_mods) if light_mods else 0,
                            len(medium_mods) if medium_mods else 0,
                            len(heavy_mods) if heavy_mods else 0))
        multiplicity_node = self.root.find('.parameterGroups/parameterGroup/'
                                           + 'multiplicity')
        multiplicity_node.text = str(multiplicity)
        max_label_node = self.root.find('.parameterGroups/parameterGroup/'
                                        + 'maxLabeledAa')
        max_label_node.text = max_label

        node = self.root.find('.parameterGroups/parameterGroup/labelMods')
        node[0].text = ';'.join(light_mods) if light_mods else ''
        if multiplicity == 3:
            MQParam._add_child(node, name='string', text=';'.join(medium_mods))
        if multiplicity > 1:
            MQParam._add_child(node, name='string',
                               text=';'.join(heavy_mods) if heavy_mods else '')

    def set_list_params(self, key, vals):
        """Set a list parameter.
        >>> t = MQParam(None, './template.xml', None)
        >>> t.set_list_params('proteases', ('test 1', 'test 2'))
        >>> len(t.root.find('.parameterGroups/parameterGroup/enzymes'))
        2
        >>> t.set_list_params('var_mods', ('Oxidation (M)', ))
        >>> t.root.find('.parameterGroups/parameterGroup/variableModifications')[0].text
        'Oxidation (M)'
        """

        params = {'var_mods' :
                  '.parameterGroups/parameterGroup/variableModifications',
                  'fixed_mods' :
                  '.parameterGroups/parameterGroup/fixedModifications',
                  'proteases' :
                  '.parameterGroups/parameterGroup/enzymes'}
        if key in params:
            node = self.root.find(params[key])
            node.clear()
            node.tag = params[key].split('/')[-1]
            for e in vals:
                MQParam._add_child(node, name='string', text=e)
        else:
            raise ValueError("Parameter {} not found.".format(key))

    def write(self):
        rough_string = ET.tostring(self.root, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        pretty = reparsed.toprettyxml(indent="\t")
        even_prettier = re.sub(r"\n\s+\n", r"\n", pretty)
        with open(self.mqpar_out, 'w') as f:
            print(even_prettier, file=f)
