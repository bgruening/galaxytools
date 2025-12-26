"""Generic utilities for linting.

Largely derived in large part from galaxy.tool_util.lint.
"""
LEVEL_ALL = "all"
LEVEL_WARN = "warn"
LEVEL_ERROR = "error"
DEFAULT_TRAINING_LINT = None


class LintContext:
    """Track running status (state) of linting."""

    def __init__(self, level=LEVEL_WARN, training_topic=DEFAULT_TRAINING_LINT):
        """Create LintContext with specified 'level' (currently unused)."""
        self.level = level
        self.training_topic = training_topic
        self.found_errors = False
        self.found_warns = False

        # self.valid_messages = []
        # self.info_messages = []
        self.warn_messages = []
        self.error_messages = []

    def __handle_message(self, message_list, message, *args, **kwds):
        if kwds or args:
            message = message.format(*args, **kwds)
        message_list.append(message)

    # def valid(self, message, *args, **kwds):
    #     self.__handle_message(self.valid_messages, message, *args, **kwds)

    # def info(self, message, *args, **kwds):
    #     self.__handle_message(self.info_messages, message, *args, **kwds)

    def error(self, message, *args, **kwds):
        """Track a linting error - a serious problem with the artifact preventing execution."""
        self.__handle_message(self.error_messages, message, *args, **kwds)

    def warn(self, message, *args, **kwds):
        """Track a linting warning - a deviation from best practices."""
        self.__handle_message(self.warn_messages, message, *args, **kwds)

    def print_messages(self):
        """Print error messages and update state at the end of linting."""
        for message in self.error_messages:
            self.found_errors = True
            print(f".. ERROR: {message}")

        if self.level != LEVEL_ERROR:
            for message in self.warn_messages:
                self.found_warns = True
                print(f".. WARNING: {message}")

        if self.level == LEVEL_ALL:
            for message in self.info_messages:
                print(f".. INFO: {message}")
            for message in self.valid_messages:
                print(f".. CHECK: {message}")
