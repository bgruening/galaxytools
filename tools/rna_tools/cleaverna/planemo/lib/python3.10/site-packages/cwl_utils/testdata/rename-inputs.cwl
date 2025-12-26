#!/usr/bin/env cwl-runner
id: InputSecondaryFileConformanceTest
class: CommandLineTool
cwlVersion: v1.2
doc: |
  Simple test to confirm the implementation of expressions returning a File within a CommandInputParameter.secondaryFile field.

  Use GREP to filter the result from ls to ensure we only get the secondary files in there.

  Related links:
  - Issue: https://github.com/common-workflow-language/cwltool/issues/1232
  - PR: https://github.com/common-workflow-language/cwltool/pull/1233
  - Discourse: https://cwl.discourse.group/t/ask-cwl-to-rename-a-secondary-file/72

inputs:
- id: inputWithSecondary
  type: File
  doc: |
    This input will with a secondary file `.accessory`. You could create these files (and its accessory) with:
    ```bash
    touch secondary_file_test.txt
    touch secondary_file_test.txt.accessory
    ```
  secondaryFiles:
    - |
      ${
        function resolveSecondary(base, secPattern) {
          if (secPattern[0] == '^') {
            var spl = base.split('.');
            var endIndex = spl.length > 1 ? spl.length - 1 : 1;
            return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
          }
          return base + secPattern;
        }
        return [{
            "class": "File",
            "location": inputs.accessory.location,
            "basename": resolveSecondary(self.basename, '^.accessory')
        }];
      }
- id: accessory
  type: File


arguments:
- "ls"
- $(inputs.inputWithSecondary.dirname)
- valueFrom: "|"
  shellQuote: false
- "grep"
- "secondary"

outputs:
- id: output_file
  type: stdout
stdout: result
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
