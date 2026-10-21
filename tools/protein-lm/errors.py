"""Exit codes shared by embed.py / predict.py / finetune.py.

  2  invalid or missing model path (preinstalled row's path does not exist,
     or an uploaded archive's extracted folder has no config.json)
  3  invalid archive (upload is not a valid tar/tgz, or it is empty)
  4  invalid FASTA (reserved; not currently raised -- parse_fasta_or_fail
     treats zero parsed records as "empty input", code 5, not a parse error)
  5  empty input (zero sequences in a FASTA file, or zero rows remaining
     after length filtering)
  6  ID mismatch between FASTA and labels table
  7  invalid label values (missing column, unparsable value, wrong class count, ...)
"""

import sys


class ToolError(Exception):
    def __init__(self, message: str, exit_code: int):
        super().__init__(message)
        self.exit_code = exit_code


def fail(message: str, exit_code: int) -> None:
    raise ToolError(message, exit_code)


def run_main(main_fn) -> None:
    """Call main_fn(); on ToolError, print to stderr and exit with its code."""
    try:
        main_fn()
    except ToolError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(exc.exit_code)
