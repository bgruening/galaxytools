__all__ = ['SilicosItError', 'WrongArgument']

class SilicosItError(Exception):
    """Base class for exceptions in Silicos-it code"""
    pass

class WrongArgument(SilicosItError):
    """
    Exception raised when argument to function is not of correct type.

    Attributes:
        function -- function in which error occurred
        msg      -- explanation of the error
    """
    def __init__(self, function, msg):
        self.function = function
        self.msg = msg
