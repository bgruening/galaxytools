"""
Module containing Exceptions that can be raised during param discovery
"""


class ArgumentParsingDiscoveryError(Exception):
    """
    Exception for general error encountered during param extraction
    """

    pass


class ArgParseImportNotFound(ArgumentParsingDiscoveryError):
    """
    Exception raised in case ArgParser was not imported
    """

    pass


class ArgParserNotUsed(ArgumentParsingDiscoveryError):
    """
    Exception raised in case no uses of argument parser were found
    """

    pass


class CouldNotFixNameError(ArgumentParsingDiscoveryError):
    pass
