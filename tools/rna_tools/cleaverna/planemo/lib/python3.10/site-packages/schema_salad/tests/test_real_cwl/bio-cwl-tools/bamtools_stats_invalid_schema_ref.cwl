#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bamtools:v2.4.1
- class: SoftwareRequirement
  packages:
    bamtools:
      specs: [ "http://identifiers.org/biotools/bamtools" ]
      version: [ "2.4.1" ]

inputs:

  bam_file:
    type:
      - File
    inputBinding:
      position: 2
      prefix: -in
    doc: |
      the input BAM file

outputs:

  log_file:
    type: File
    outputBinding:
      glob: "stats.log"

  total_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var totalReads = parseInt(s.substring ( s.indexOf("Totalreads")+11, s.indexOf("\t", (s.indexOf("Totalreads")))  ));
          return totalReads;
        }

  mapped_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var mappedreads = parseInt(s.substring ( s.indexOf("Mappedreads")+12, s.indexOf("\t", (s.indexOf("Mappedreads")))  ));
          return mappedreads;
        }

  forward_strand_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var forwardstrand = parseInt(s.substring ( s.indexOf("Forwardstrand")+14, s.indexOf("\t", (s.indexOf("Forwardstrand")))  ));
          return forwardstrand;
        }

  reverse_strand_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var reversestrand = parseInt(s.substring ( s.indexOf("Reversestrand")+14, s.indexOf("\t", (s.indexOf("Reversestrand")))  ));
          return reversestrand;
        }

  failed_QC_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var failedQC = parseInt(s.substring ( s.indexOf("FailedQC")+9, s.indexOf("\t", (s.indexOf("FailedQC")))  ));
          return failedQC;
        }

  duplicate_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var duplicates = parseInt(s.substring ( s.indexOf("Duplicates")+11, s.indexOf("\t", (s.indexOf("Duplicates")))  ));
          return duplicates;
        }

  pairedend_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: "stats.log"
      outputEval: |
        ${
          var s = self[0].contents.replace(/ /g,'').replace(/ *\([^)]*\) */g,'');
          var pairedendreads = parseInt(s.substring ( s.indexOf("Paired-endreads")+16, s.indexOf("\t", (s.indexOf("Paired-endreads")))  ));
          return pairedendreads;
        }

baseCommand: [bamtools, stats]
arguments:
  - valueFrom: $('> ' + 'stats.log')
    position: 1000
    shellQuote: false

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/9.0/schemaorg-current-http.rdf

s:name: "bamtools_stats"
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Tool runs `bamtools stats' to calculate general alignment statistics from the input BAM file

  `-insert` parameter is not implemented