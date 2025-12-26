from pathlib import Path
from typing import (
    Optional,
    Union,
)

import yaml
from pydantic import (
    BaseModel,
    Extra,
    RootModel,
)

StrOrPath = Union[Path, str]


class RepositoryInstallTarget(BaseModel):
    name: str
    owner: str
    tool_shed_url: Optional[str] = None
    tool_panel_section_id: Optional[str] = None
    tool_panel_section_label: Optional[str] = None
    revisions: Optional[list[str]] = None
    install_tool_dependencies: Optional[bool] = None
    install_repository_dependencies: Optional[bool] = None
    install_resolver_dependencies: Optional[bool] = None


class RepositoryInstallTargets(BaseModel):
    """ """

    api_key: Optional[str] = None
    galaxy_instance: Optional[str] = None
    tools: list[RepositoryInstallTarget]


class DataManager(BaseModel, extra=Extra.forbid):
    tags: list[str]
    tool_id: str


class DataManagers(RootModel):
    root: dict[str, DataManager]


class Genome(BaseModel):
    id: str  # The unique id of the data in Galaxy
    description: Optional[str] = None  # The description of the data, including its taxonomy, version and date
    dbkey: Optional[str] = None
    source: Optional[str] = (
        None  # The source of the data. Can be: 'ucsc', an NCBI accession number or a URL to a fasta file.
    )

    # The following fields are currently purely for human consumption and unused by
    # IDC infrastructure.
    doi: Optional[str] = None  # Any DOI associated with the data
    blob: Optional[str] = None  # A blob for any other pertinent information
    checksum: Optional[str] = None  # A SHA256 checksum of the original
    version: Optional[str] = None  # Any version information associated with the data

    # Description of actions (data managers) to run on target genome.
    indexers: Optional[
        list[str]
    ]  # indexers to run - keyed on repository name - see data_managers.yml for how to resolve these to tools
    skiplist: Optional[list[str]] = (
        None  # unimplemented: but if we implement classes of indexers, these will be ones to skip
    )


class Genomes(BaseModel):
    genomes: list[Genome]


def _read_yaml(path: StrOrPath):
    with open(path) as f:
        return yaml.safe_load(f)


def read_data_managers(path: StrOrPath) -> DataManagers:
    return DataManagers(root=_read_yaml(path))


def read_genomes(path: StrOrPath) -> Genomes:
    return Genomes(**_read_yaml(path))


def read_tools(path: StrOrPath) -> RepositoryInstallTargets:
    return RepositoryInstallTargets(**_read_yaml(path))
