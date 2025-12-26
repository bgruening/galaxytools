#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  tarfile:
    type: File
    inputBinding:
      prefix: --file
  extractfile:
    type: string
    inputBinding:
      position: 1
baseCommand: [tar, --extract]
outputs:
  extracted_file:
    type: File
    outputBinding:
      glob: $(inputs.extractfile)
