"""This file contains the parser for shed_tools"""

import argparse

from .common_parser import (
    get_common_args,
    HideUnderscoresHelpFormatter,
)


def parser():
    """construct the parser object"""
    common_arguments = get_common_args(log_file=True)
    shed_parser = argparse.ArgumentParser()
    subparsers = shed_parser.add_subparsers()

    # A list of defaults is needed. Otherwise the shed-tools install parser will not return
    # update_tools in the name space and shed-tool update will not return all the install
    # variables.
    shed_parser.set_defaults(
        action="install",
        tool_list_file=None,
        tool_yaml=None,
        owner=None,
        name=None,
        tool_panel_section_id=None,
        tool_panel_section_label=None,
        revisions=None,
        tool_shed_url=None,
        skip_tool_dependencies=False,
        install_resolver_dependencies=False,
        force_latest_revision=False,
        test=False,
        test_user_api_key=None,
        test_user="ephemeris@galaxyproject.org",
        test_json="tool_test_output.json",
        test_existing=False,
        parallel_tests=1,
        client_test_config=None,
    )

    # SUBPARSERS
    install_command_parser = subparsers.add_parser(
        "install",
        help="This installs tools in Galaxy from the Tool Shed." "Use shed-tools install --help for more information",
        formatter_class=HideUnderscoresHelpFormatter,
        parents=[common_arguments],
    )
    update_command_parser = subparsers.add_parser(
        "update",
        help="This updates all tools in Galaxy to the latest revision. "
        "Use shed-tools update --help for more information",
        formatter_class=HideUnderscoresHelpFormatter,
        parents=[common_arguments],
    )

    test_command_parser = subparsers.add_parser(
        "test",
        help="This tests the supplied list of tools in Galaxy. " "Use shed-tools test --help for more information",
        formatter_class=HideUnderscoresHelpFormatter,
        parents=[common_arguments],
    )

    # SUBPARSER DEFAULTS
    update_command_parser.set_defaults(action="update")

    test_command_parser.set_defaults(action="test")
    install_command_parser.set_defaults(action="install")

    # COMMON OPTIONS
    for command_parser in [
        update_command_parser,
        install_command_parser,
        test_command_parser,
    ]:
        command_parser.add_argument(
            "-t",
            "--tools-file",
            "--toolsfile",
            dest="tool_list_file",
            help="Tools file to use (see tool_list.yaml.sample)",
        )
        command_parser.add_argument(
            "-y",
            "--yaml-tool",
            "--yaml_tool",
            dest="tool_yaml",
            help="Install tool represented by yaml string",
        )
        command_parser.add_argument(
            "--name",
            help="The name of the tool to install (only applicable " "if the tools file is not provided).",
        )
        command_parser.add_argument(
            "--owner",
            help="The owner of the tool to install (only applicable " "if the tools file is not provided).",
        )
        command_parser.add_argument(
            "--revisions",
            default=None,
            nargs="*",
            dest="revisions",
            help="The revisions of the tool repository that will be installed. "
            "All revisions must be specified after this flag by a space."
            "Example: --revisions 0a5c7992b1ac f048033da666"
            "(Only applicable if the tools file is not provided).",
        )
        command_parser.add_argument(
            "--tool-shed",
            "--toolshed",
            dest="tool_shed_url",
            help="The Tool Shed URL where to install the tool from. "
            "This is applicable only if the tool info is "
            "provided as an option vs. in the tools file.",
        )

    # OPTIONS COMMON FOR UPDATE AND INSTALL

    for command_parser in [update_command_parser, install_command_parser]:
        command_parser.add_argument(
            "--install-tool-dependencies",
            "--install_tool_dependencies",
            action="store_true",
            dest="install_tool_dependencies",
            default=False,
            help="Turn on installation of tool dependencies using classic toolshed packages. "
            "Can be overwritten on a per-tool basis in the tools file.",
        )
        command_parser.add_argument(
            "--skip-install-resolver-dependencies",
            "--skip_install_resolver_dependencies",
            action="store_false",
            dest="install_resolver_dependencies",
            default=False,
            help="Skip installing tool dependencies through resolver (e.g. conda). "
            "Will be ignored on galaxy releases older than 16.07. "
            "Can be overwritten on a per-tool basis in the tools file",
        )
        command_parser.add_argument(
            "--skip-install-repository-dependencies",
            "--skip_install_repository_dependencies",
            action="store_false",
            dest="install_repository_dependencies",
            default=False,
            help="Skip installing the repository dependencies.",
        )
        command_parser.add_argument(
            "--test",
            action="store_true",
            dest="test",
            help="Run tool tests on install tools, requires Galaxy 18.05 or newer.",
        )
        command_parser.add_argument(
            "--test-existing",
            "--test_existing",
            action="store_true",
            help="If testing tools during install, also run tool tests on repositories already installed "
            "(i.e. skipped repositories).",
        )
        command_parser.add_argument(
            "--test-json",
            "--test_json",
            dest="test_json",
            default="tool_test_output.json",
            help="If testing tools, record tool test output to specified file. "
            "This file can be turned into reports with ``planemo test_reports <output.json>``.",
        )
        command_parser.add_argument(
            "--test-user-api-key",
            "--test_user_api_key",
            dest="test_user",
            help="If testing tools, a user is needed to execute the tests. "
            "This can be different the --api_key which is assumed to be an admin key. "
            "If --api_key is a valid user (e.g. it is not a master API key) this does "
            "not need to be specified and --api_key will be reused.",
        )
        command_parser.add_argument(
            "--test-user",
            "--test_user",
            dest="test_user",
            help="If testing tools, a user is needed to execute the tests. "
            "If --api_key is a master api key (i.e. not tied to a real user) and "
            "--test_user_api_key isn't specified, this user email will be used. This "
            "user will be created if needed.",
        )
        command_parser.add_argument(
            "--parallel-tests",
            "--parallel_tests",
            dest="parallel_tests",
            default=1,
            type=int,
            help="Specify the maximum number of tests that will be run in parallel.",
        )

    # OPTIONS UNIQUE TO INSTALL

    install_command_parser.add_argument(
        "--section",
        dest="tool_panel_section_id",
        help="Galaxy tool panel section ID where the tool will "
        "be installed (the section must exist in Galaxy; "
        "only applicable if the tools file is not provided).",
    )
    install_command_parser.add_argument(
        "--section-label",
        "--section_label",
        default=None,
        dest="tool_panel_section_label",
        help="Galaxy tool panel section label where tool will be installed "
        "(if the section does not exist, it will be created; "
        "only applicable if the tools file is not provided).",
    )

    install_command_parser.add_argument(
        "--latest",
        action="store_true",
        dest="force_latest_revision",
        help="Will override the revisions in the tools file and always install the latest revision.",
    )

    # OPTIONS UNIQUE TO TEST
    # Same test_json as above but language modified for test instead of install/update.
    test_command_parser.add_argument(
        "--test-json",
        "--test_json",
        default="tool_test_output.json",
        dest="test_json",
        help="Record tool test output to specified file. "
        "This file can be turned into reports with ``planemo test_reports <output.json>``.",
    )

    test_command_parser.add_argument(
        "--test-user-api-key",
        "--test_user_api_key",
        dest="test_user_api_key",
        help="A user is needed to execute the tests. "
        "This can be different the --api_key which is assumed to be an admin key. "
        "If --api_key is a valid user (e.g. it is not a master API key) this does "
        "not need to be specified and --api_key will be reused.",
    )
    test_command_parser.add_argument(
        "--test-user",
        "--test_user",
        dest="test_user",
        help="A user is needed to execute the tests. "
        "If --api_key is a master api key (i.e. not tied to a real user) and "
        "--test_user_api_key isn't specified, this user email will be used. This "
        "user will be created if needed.",
    )
    test_command_parser.add_argument(
        "--test-history-name",
        "--test_history_name",
        dest="test_history_name",
        default=None,
        help="Use existing history or create history with provided name if none exists. "
        "If --test_history_name is not set, a new history with a default name will always "
        "be created. If multiple histories match the provided name, the first (newest) "
        "one returned by the Galaxy API will be selected.",
    )
    test_command_parser.add_argument(
        "--parallel-tests",
        "--parallel_tests",
        dest="parallel_tests",
        default=1,
        type=int,
        help="Specify the maximum number of tests that will be run in parallel.",
    )
    test_command_parser.add_argument(
        "--test-all-versions",
        "--test_all_versions",
        action="store_true",
        dest="test_all_versions",
        help="Run tests on all installed versions of tools.  This will only "
        "apply for tools where revisions have not been provided through "
        "the --revisions arg, --tool_file or --tool_yaml.",
    )
    test_command_parser.add_argument(
        "--client-test-config",
        "--client_test_config",
        dest="client_test_config",
        help="Annotate expectations about tools in client testing YAML " "configuration file.",
    )

    return shed_parser
