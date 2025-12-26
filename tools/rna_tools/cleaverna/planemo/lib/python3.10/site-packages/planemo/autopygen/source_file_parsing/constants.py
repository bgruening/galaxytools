from sys import version_info

# stdlib_module_names is only available from python 3.10
if version_info >= (3, 10):
    from sys import stdlib_module_names  # type: ignore[attr-defined]
else:
    from stdlib_list import stdlib_list

    stdlib_module_names = stdlib_list()

WARNING_STRING = "##!_FIXME_!##"

STD_LIB_MODULE_NAMES = frozenset(stdlib_module_names)
