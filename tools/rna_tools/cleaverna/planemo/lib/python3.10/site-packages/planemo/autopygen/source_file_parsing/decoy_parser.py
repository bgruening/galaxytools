import argparse
import ast
from typing import Optional


class DecoyParser(argparse.ArgumentParser):
    """
    Decoy class that is injected into code that initializes parser
    It removes dependency on the internal representation of actions
    """

    class Subparsers:
        def __init__(self):
            self.parsers = []

        def add_parser(self, name, **kwargs):
            dp = DecoyParser(name)
            self.parsers.append(dp)
            return dp

    class Action:
        def __init__(self, section, argument, action, kwargs):
            self.argument = argument

            if action is None:
                action = "STORE"

            self.action = action.upper()
            self.kwargs = kwargs
            self.scope = section

    class Section:
        def __init__(self, parent=None, name=None, description=None):
            self.name = name
            self.description = description
            self.parent = parent
            self.actions = []
            self.subsections = []

        def get_actions_recursive(self):
            for action in self.actions:
                yield action

            for subsection in self.subsections:
                yield from subsection.get_actions_recursive()

    def __init__(self, name="default"):
        self.name = name
        self.default_section = self.Section(name=name)
        self.sub_parsers = []
        super().__init__()

    def report_arguments_and_groups(self):
        yield from self.default_section.get_actions_recursive()

    def save_action(self, section, *args, **kwargs):
        section.actions.append(self.Action(section, args[-1], kwargs.get("action", "STORE"), kwargs))

    def add_argument(self, *args, **kwargs):
        self.save_action(self.default_section, *args, **kwargs)
        return super().add_argument(*args, **kwargs)

    def add_argument_for_arg_group(self, section, group):
        def add_argument(*args, **kwargs):
            self.save_action(section, *args, **kwargs)
            return super(type(group), group).add_argument(*args, **kwargs)

        return add_argument

    def create_custom_add_argument_group(self, original, parent_section):
        def custom_add_argument_group(*args, **kwargs):
            from copy import copy

            new_grp = original.add_argument_group(*args, **kwargs)
            if len(args) == 0:
                subsection = self.Section(parent_section, **kwargs)
            else:
                subsection = self.Section(parent_section, *args)
            parent_section.subsections.append(subsection)

            new_grp.add_argument = self.add_argument_for_arg_group(subsection, new_grp)
            new_grp.add_argument_group = self.create_custom_add_argument_group(copy(new_grp), subsection)
            return new_grp

        return custom_add_argument_group

    def add_argument_group(self, *args, **kwargs):
        from copy import copy

        arg_group = super().add_argument_group(*args, **kwargs)

        if len(args) == 0:
            subsection = self.Section(self.default_section, **kwargs)
        else:
            subsection = self.Section(self.default_section, *args)

        arg_group.add_argument = self.add_argument_for_arg_group(subsection, arg_group)
        arg_group.add_argument_group = self.create_custom_add_argument_group(copy(arg_group), subsection)
        self.default_section.subsections.append(subsection)
        return arg_group

    def add_subparsers(self, **kwargs):
        mock_subparser = self.Subparsers()
        self.sub_parsers.append(mock_subparser)
        return mock_subparser


def obtain_class_def() -> Optional[ast.ClassDef]:
    file = open(__file__, "r")
    module = ast.parse(file.read())
    file.close()
    return next((item for item in module.body if type(item) is ast.ClassDef and item.name == "DecoyParser"), None)


class CustomParserUnavailableException(Exception):
    pass
