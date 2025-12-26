"""YAML loading utilities for gxformat2."""
from collections import OrderedDict

try:
    from galaxy.model.custom_types import MutationDict  # type: ignore
except ImportError:
    MutationDict = None
import yaml


def ordered_load_path(path: str, **kwds):
    """Safe and ordered load of YAML from specified path."""
    with open(path) as f:
        return ordered_load(f, **kwds)


def ordered_load(stream, Loader=yaml.SafeLoader, **kwds):
    """Safe and ordered load of YAML from stream."""
    class OrderedLoader(Loader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return OrderedDict(loader.construct_pairs(node))

    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)

    return yaml.load(stream, OrderedLoader, **kwds)


def ordered_dump(data, stream=None, Dumper=yaml.SafeDumper, **kwds):
    """Safe and ordered dump of YAML to stream."""
    class OrderedDumper(Dumper):
        pass

    def _dict_representer(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
            list(data.items()))
    OrderedDumper.add_representer(OrderedDict, _dict_representer)
    if MutationDict is not None:
        OrderedDumper.add_representer(MutationDict, _dict_representer)

    return yaml.dump(data, stream, OrderedDumper, **kwds)


def ordered_dump_to_path(as_dict: dict, path: str):
    """Safe and ordered dump of YAML to path."""
    with open(path, "w") as f:
        ordered_dump(as_dict, f)


__all__ = (
    'ordered_dump',
    'ordered_dump_to_path',
    'ordered_load',
    'ordered_load_path',
)
