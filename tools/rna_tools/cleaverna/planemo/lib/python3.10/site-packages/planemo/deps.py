"""Abstractions for building dependency resolution configurations."""

import tempfile
from string import Template

import click

from planemo.conda import build_conda_context

CONDA_DEPENDENCY_RESOLUTION_CONF = """<dependency_resolvers>
  <conda ${attributes} />
  <conda versionless="true" ${attributes} />
</dependency_resolvers>
"""

# Like Conda resolution above, but allow tool shed packages to be used for
# shed_serve and shed_test.
DEFAULT_DEPENDENCY_RESOLUTION_CONF = """<dependency_resolvers>
  <tool_shed_packages />
  <conda ${attributes} />
  <conda versionless="true" ${attributes} />
</dependency_resolvers>
"""

NO_DEPENDENCY_RESOLUTION_CONF = """<dependency_resolvers>
</dependency_resolvers>
"""

BREW_DEPENDENCY_RESOLUTION_CONF = """<dependency_resolvers>
  <homebrew />
  <!--
  <homebrew versionless="true" />
  -->
</dependency_resolvers>
"""

SHED_DEPENDENCY_RESOLUTION_CONF = """<dependency_resolvers>
  <tool_shed_tap />
</dependency_resolvers>
"""

# Provide some shortcuts for simple/common dependency resolutions strategies.
STOCK_DEPENDENCY_RESOLUTION_STRATEGIES = {
    "brew_dependency_resolution": BREW_DEPENDENCY_RESOLUTION_CONF,
    "shed_dependency_resolution": SHED_DEPENDENCY_RESOLUTION_CONF,
    "conda_dependency_resolution": CONDA_DEPENDENCY_RESOLUTION_CONF,
    "no_dependency_resolution": NO_DEPENDENCY_RESOLUTION_CONF,
    "default_dependency_resolution": DEFAULT_DEPENDENCY_RESOLUTION_CONF,
}


def ensure_dependency_resolvers_conf_configured(ctx, kwds, resolvers_conf=None):
    """Use supplied CLI options (kwds) to find or build a dependency resolvers file.

    Set new path in kwds if needed.
    """
    _validate_dependency_resolution_options(kwds)
    always_specify_attribute = object()

    dependency_attribute_kwds = {
        "conda_prefix": None,
        "conda_exec": None,
        "conda_debug": False,
        "conda_copy_dependencies": False,
        "conda_auto_init": always_specify_attribute,
        "conda_auto_install": always_specify_attribute,
        "conda_ensure_channels": "",
        "conda_use_local": False,
    }
    attributes = []

    def add_attribute(key, value):
        attributes.append(f'{key}="{value}"')

    conda_prefix_specified = False
    for key, default_value in dependency_attribute_kwds.items():
        value = kwds.get(key, default_value)
        if value != default_value:
            conda_prefix_specified = conda_prefix_specified or (key == "conda_prefix")
            # Strip leading prefix (conda_) off attributes
            attribute_key = "_".join(key.split("_")[1:])
            add_attribute(attribute_key, value)

    conda_context = build_conda_context(ctx, **kwds)
    if not conda_prefix_specified:
        add_attribute("prefix", conda_context.conda_prefix)
    add_attribute("condarc_override", conda_context.condarc_override)

    attribute_str = " ".join(attributes)

    if kwds.get("dependency_resolvers_config_file", None):
        resolution_type = "__explicit__"
    else:
        resolution_type = "default_dependency_resolution"
        for key in STOCK_DEPENDENCY_RESOLUTION_STRATEGIES:
            if kwds.get(key):
                resolution_type = key

    if resolution_type != "__explicit__":
        # Planemo manages the dependency resolve conf file.
        template_str = STOCK_DEPENDENCY_RESOLUTION_STRATEGIES[resolution_type]
        conf_contents = Template(template_str).safe_substitute({"attributes": attribute_str})
        if resolvers_conf is None:
            resolvers_conf = tempfile.NamedTemporaryFile(delete=False).name
        with open(resolvers_conf, "w") as fh:
            fh.write(conf_contents)
        ctx.vlog(
            "Writing dependency_resolvers_config_file to path %s with contents [%s]",
            resolvers_conf,
            conf_contents,
        )
        kwds["dependency_resolvers_config_file"] = resolvers_conf


def _validate_dependency_resolution_options(kwds):
    resolutions_strategies = [
        "brew_dependency_resolution",
        "dependency_resolvers_config_file",
        "shed_dependency_resolution",
        "conda_dependency_resolution",
    ]

    selected_strategies = 0
    for key in resolutions_strategies:
        if kwds.get(key):
            selected_strategies += 1

    if selected_strategies > 1:
        message = "At most one option from [%s] may be specified"
        raise click.UsageError(message % resolutions_strategies)
