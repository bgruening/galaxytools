"""
Example from README.md.

Please synchronize all changes between the two.
"""

# SPDX-License-Identifier: Apache-2.0
from pathlib import Path

from cwl_utils.parser import load_document_by_uri, save

# File Input - This is the only thing you will need to adjust or take in as an input to your function:
cwl_file = Path("testdata/md5sum.cwl")  # or a plain string works as well

# Import CWL Object
cwl_obj = load_document_by_uri(cwl_file)

# View CWL Object
print("List of object attributes:\n{}".format("\n".join(map(str, dir(cwl_obj)))))

# Export CWL Object into a built-in typed object
saved_obj = save(cwl_obj)
print(f"Export of the loaded CWL object: {saved_obj}.")
