class: Workflow
cwlVersion: draft-3
inputs:
- id: '#input_file'
  type: ['null', File]
outputs:
- id: '#validatefiles_report'
  source: '#validatefiles.report'
  type: ['null', File]
- id: '#md5_report'
  source: '#md5.report'
  type: ['null', File]
requirements:
- {class: InlineJavascriptRequirement}
steps:
- id: '#md5'
  inputs:
  - {id: '#md5.input_file', source: '#input_file'}
  outputs:
  - {id: '#md5.report'}
  run: md5.cwl
- id: '#validatefiles'
  inputs:
  - {id: '#validatefiles.input_file', source: '#input_file'}
  - {id: '#validatefiles.type'}
  outputs:
  - {id: '#validatefiles.report'}
  run: validate.cwl
