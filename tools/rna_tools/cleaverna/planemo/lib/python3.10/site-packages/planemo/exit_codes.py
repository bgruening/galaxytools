"""Numeric constants for planemo exit codes."""

# Operation succeeded.
EXIT_CODE_OK = 0

# Generic failure.
EXIT_CODE_GENERIC_FAILURE = 1

# If iterating or tools for instance, not tool xml files are found.
EXIT_CODE_NO_SUCH_TARGET = 2

# An unknown file type was attempted to be processed with planemo.
EXIT_CODE_UNKNOWN_FILE_TYPE = 4

# An unsupported file type was supplied for a given operation.
EXIT_CODE_UNSUPPORTED_FILE_TYPE = 5

# A dependency of this operation was unavailable (e.g. conda).
EXIT_CODE_FAILED_DEPENDENCIES = 6

# Attempt to do a one time action that already happened.
EXIT_CODE_ALREADY_EXISTS = 7


class ExitCodeException(Exception):
    """Exception used by planemo framework to track exit codes for CLI."""

    def __init__(self, exit_code):
        """Specify integer exit code to exit planemo with."""
        self.exit_code = exit_code
