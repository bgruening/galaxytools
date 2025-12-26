"""Module describing the planemo ``conda_env`` command."""

import click
from galaxy.tool_util.deps import conda_util

from planemo import options
from planemo.cli import command_function
from planemo.conda import (
    build_conda_context,
    collect_conda_targets,
)
from planemo.io import (
    error,
    ps1_for_path,
)

SOURCE_COMMAND = """
PRE_CONDA_PS1=$PS1
source %s %s
if [[ -n $BASH_VERSION ]]; then
    hash -r
elif [[ -n $ZSH_VERSION ]]; then
    rehash
else
echo 'Only bash and zsh are supported'
    return 1
fi
PS1="%s"
echo 'Deactivate environment with conda_env_deactivate'
alias conda_env_deactivate="source %s; %s"
"""


@click.command("conda_env")
@options.optional_tools_arg()
@options.conda_target_options()
# @options.skip_install_option()  # TODO
@command_function
def cli(ctx, path, **kwds):
    """Activate a conda environment for tool.

    Source the output of this command to activate a conda environment for this
    tool.

    \b
        $ . <(planemo conda_env seqtk_seq.xml)
        Deactivate environment with conda_env_deactivate
        (seqtk_seq_v6) $ which seqtk
        /home/planemo/miniconda2/envs/jobdepsDkzcjjfecc6d406196737781ff4456ec60975c137e04884e4f4b05dc68192f7cec4656/bin/seqtk
        (seqtk_seq_v6) $ conda_env_deactivate
        $

    """
    conda_context = build_conda_context(ctx, use_planemo_shell_exec=False, **kwds)
    conda_targets = collect_conda_targets(ctx, [path])
    installed_conda_targets = conda_util.filter_installed_targets(conda_targets, conda_context=conda_context)
    env_name, exit_code = conda_util.build_isolated_environment(
        installed_conda_targets, conda_context=conda_context, quiet=True
    )
    if exit_code:
        error("Failed to build environment for request.")
        return 1

    ps1 = ps1_for_path(path, base="PRE_CONDA_PS1")
    remove_env = f"{conda_context.conda_exec} env remove -y --name '{env_name}'"
    deactivate = conda_context.deactivate
    activate = conda_context.activate
    command = SOURCE_COMMAND % (activate, env_name, ps1, deactivate, remove_env)
    print(command)
