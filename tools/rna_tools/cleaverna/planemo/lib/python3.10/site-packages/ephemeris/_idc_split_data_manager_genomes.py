#!/usr/bin/env python
"""Helper script for IDC - not yet meant for public consumption.

This script splits genomes.yml into tasks that are meant to be sent to
run_data_managers.py - while excluding data managers executions specified
by genomes.yml that have already been executed and appear in the target
installed data table configuration.
"""
import logging
import os
import re
import xml.etree.ElementTree as ElementTree
from copy import deepcopy
from typing import (
    Any,
    Callable,
    Optional,
)

import requests
import yaml
from galaxy.util import safe_makedirs
from pydantic import (
    BaseModel,
    Extra,
    RootModel,
)

from . import get_galaxy_connection
from ._idc_data_managers_to_tools import (
    DataManager,
    read_data_managers_configuration,
)
from .common_parser import get_common_args
from .ephemeris_log import (
    disable_external_library_logging,
    setup_global_logger,
)

IsBuildComplete = Callable[[str, str], bool]
TASK_FILE_NAME = "run_data_managers.yaml"
DEFAULT_TOOL_ID_MODE = "tool_shed_guid"
UCSC_DSN_URL = "http://genome.cse.ucsc.edu/cgi-bin/das/dsn"

log = logging.getLogger(__name__)


class Filters:
    stage: Optional[int] = None
    data_manager: Optional[str] = None
    build_id: Optional[str] = None

    def filter_out_data_manager(self, data_manager: str) -> bool:
        return bool(self.data_manager and data_manager != self.data_manager)

    def filter_out_build_id(self, build_id: str) -> bool:
        return bool(self.build_id and build_id != self.build_id)

    def filter_out_stage(self, stage: int) -> bool:
        return bool(self.stage is not None and self.stage != stage)


class SplitOptions:
    merged_genomes_path: str
    split_genomes_path: str
    data_managers_path: str
    is_build_complete: IsBuildComplete
    tool_id_mode: str = DEFAULT_TOOL_ID_MODE
    filters: Filters = Filters()


def tool_id_for(indexer: str, data_managers: dict[str, DataManager], mode: str) -> str:
    data_manager = data_managers[indexer]
    assert data_manager, f"Could not find a target data manager for indexer name {indexer}"
    tool_shed_guid = data_manager.tool_id
    if mode == "short":
        _ts, _, _owner, _repo_name, rest = tool_shed_guid.split("/", 4)
        if "/" in rest:
            print(rest)
            return rest.split("/")[0]
        else:
            return rest
    else:
        return tool_shed_guid


class RunDataManager(BaseModel):
    id: str
    items: Optional[list[Any]] = None
    params: Optional[list[Any]] = None
    data_table_reload: Optional[list[str]] = None


class RunDataManagers(BaseModel):
    data_managers: list[RunDataManager]


class DataManager(BaseModel, extra=Extra.forbid):
    tags: list[str]
    tool_id: str


class DataManagers(RootModel):
    root: dict[str, DataManager]


class Genome(BaseModel):
    pass


class Genomes(BaseModel):
    genomes: list[Genome]


def ucsc_description_for_build(requested_build: str) -> str:
    # from galaxy/cron/parse_builds.py
    url = UCSC_DSN_URL
    text = requests.get(url).text
    tree = ElementTree.fromstring(text)

    for dsn in tree:
        build = dsn.find("SOURCE").attrib["id"]
        if build != requested_build:
            continue

        description = dsn.find("DESCRIPTION").text.replace(" - Genome at UCSC", "").replace(" Genome at UCSC", "")

        fields = description.split(" ")
        temp = fields[0]
        for i in range(len(fields) - 1):
            if temp == fields[i + 1]:
                fields.pop(i + 1)
            else:
                temp = fields[i + 1]
        description = " ".join(fields)
        return description

    raise Exception(f"Could not fetch UCSC description for build, a description must be provided: {requested_build}")


def write_run_data_manager_to_file(run_data_manager: RunDataManager, path: str):
    parent, _ = os.path.split(path)
    if not os.path.exists(parent):
        safe_makedirs(parent)
    run_data_managers = RunDataManagers(data_managers=[run_data_manager])
    with open(path, "w") as of:
        yaml.safe_dump(run_data_managers.dict(exclude_unset=True), of)


def walk_over_incomplete_runs(split_options: SplitOptions):
    data_managers = read_data_managers_configuration(split_options.data_managers_path)
    with open(split_options.merged_genomes_path) as f:
        genomes_all = yaml.safe_load(f)
    genomes = genomes_all["genomes"]
    for genome in genomes:
        build_id = genome["id"]
        if split_options.filters.filter_out_build_id(build_id):
            continue

        fetch_indexer = "data_manager_fetch_genome_dbkeys_all_fasta"
        do_fetch = not split_options.filters.filter_out_data_manager(fetch_indexer)
        source = genome.get("source")
        if source is None:
            do_fetch = False
        if do_fetch and split_options.filters.filter_out_stage(0):
            do_fetch = False

        if do_fetch and not split_options.is_build_complete(build_id, fetch_indexer):
            log.info(f"Fetching: {build_id}")
            fetch_tool_id = tool_id_for(fetch_indexer, data_managers, split_options.tool_id_mode)
            fetch_params = []
            fetch_params.append({"dbkey_source|dbkey_source_selector": "new"})
            fetch_params.append({"dbkey_source|dbkey": genome["id"]})
            description = genome.get("description")
            source = genome.get("source")
            if source == "ucsc":
                if not description:
                    description = ucsc_description_for_build(genome["id"])
                fetch_params.append({"reference_source|reference_source_selector": "ucsc"})
                fetch_params.append({"reference_source|requested_dbkey": genome["id"]})
                fetch_params.append({"sequence_name": description})
            elif re.match("^[A-Z_]+[0-9.]+", source):
                fetch_params.append({"reference_source|reference_source_selector": "ncbi"})
                fetch_params.append({"reference_source|requested_identifier": source})
                fetch_params.append({"sequence_name": genome["description"]})
                fetch_params.append({"sequence.id": genome["id"]})
            elif re.match("^http", source):
                fetch_params.append({"reference_source|reference_source_selector": "url"})
                fetch_params.append({"reference_source|user_url": source})
                fetch_params.append({"sequence_name": genome["description"]})
                fetch_params.append({"sequence.id": genome["id"]})

            if description:
                fetch_params.append({"dbkey_source|dbkey_name": description})

            fetch_run_data_manager = RunDataManager(
                id=fetch_tool_id,
                params=fetch_params,
                # Not needed according to Marius
                # data_table_reload=["all_fasta", "__dbkeys__"],
            )
            yield (build_id, fetch_indexer, fetch_run_data_manager)
        else:
            log.debug(f"Fetch is already completed: {build_id}")

        indexers = genome.get("indexers", [])
        for indexer in indexers:
            if split_options.filters.filter_out_data_manager(indexer):
                continue

            if split_options.filters.filter_out_stage(1):
                continue

            if split_options.is_build_complete(build_id, indexer):
                log.debug(f"Build is already completed: {build_id} {indexer}")
                continue

            log.info(f"Building: {build_id} {indexer}")

            tool_id = tool_id_for(indexer, data_managers, split_options.tool_id_mode)
            params = [
                {"all_fasta_source": "{{ item.id }}"},
                {"sequence_name": "{{ item.name }}"},
                {"sequence_id": "{{ item.id }}"},
            ]
            # why is this not pulled from the data managers conf? -nate
            if re.search("bwa", tool_id):
                params.append({"index_algorithm": "bwtsw"})
            if re.search("color_space", tool_id):
                continue

            item = deepcopy(genome)
            item.pop("indexers", None)
            item.pop("skiplist", None)

            run_data_manager = RunDataManager(
                id=tool_id,
                params=params,
                items=[item],
            )
            yield (build_id, indexer, run_data_manager)


def split_genomes(split_options: SplitOptions) -> None:
    def write_task_file(build_id: str, indexer: str, run_data_manager: RunDataManager):
        split_genomes_path = split_options.split_genomes_path
        if not os.path.exists(split_options.split_genomes_path):
            safe_makedirs(split_genomes_path)

        task_file_dir = os.path.join(split_genomes_path, build_id, indexer)
        task_file = os.path.join(task_file_dir, TASK_FILE_NAME)
        write_run_data_manager_to_file(run_data_manager, task_file)

    for build_id, indexer, run_data_manager in walk_over_incomplete_runs(split_options):
        write_task_file(build_id, indexer, run_data_manager)


class GalaxyHistoryIsBuildComplete:
    def __init__(self, history_names: list[str]):
        self._history_names = history_names

    def __call__(self, build_id: str, indexer_name: str) -> bool:
        target_history_name = f"idc-{build_id}-{indexer_name}"
        return target_history_name in self._history_names


class CVMFSPublishIsComplete:
    def __init__(self, records: dict[str, list[str]]):
        self.records = records

    def __call__(self, build_id: str, indexer_name: str) -> bool:
        return indexer_name in self.records.get(build_id, [])


def _parser():
    """returns the parser object."""
    # login required to check history...
    parser = get_common_args(login_required=True, log_file=True)
    parser.add_argument("--merged-genomes-path", "-m", default="genomes.yml")
    parser.add_argument("--split-genomes-path", "-s", default="data_manager_tasks")
    parser.add_argument("--data-managers-path", default="data_managers.yml")
    parser.add_argument("--complete-check-cvmfs", default=False, action="store_true")
    parser.add_argument("--cvmfs-root", default="/cvmfs/idc.galaxyproject.org")

    parser.add_argument("--tool-id-mode", choices=["tool_shed_guid", "short"], default=DEFAULT_TOOL_ID_MODE)

    # filters
    parser.add_argument("--filter-stage", default=None)
    parser.add_argument("--filter-data-manager", default=None)
    parser.add_argument("--filter-build-id", default=None)

    return parser


def get_galaxy_history_names(args) -> list[str]:
    gi = get_galaxy_connection(args, login_required=True)
    return [h["name"] for h in gi.histories.get_histories()]


def get_regular_files(dirname: str) -> list[str]:
    return [f for f in os.listdir(dirname) if not f.startswith(".")]


def get_cvmfs_publish_records(args) -> dict[str, list[str]]:
    records = {}
    records_dir = os.path.join(args.cvmfs_root, "record")
    for build_id in get_regular_files(records_dir):
        records[build_id] = get_regular_files(os.path.join(records_dir, build_id))
    return records


def main():
    disable_external_library_logging()
    parser = _parser()
    args = parser.parse_args()
    log = setup_global_logger(name=__name__, log_file=args.log_file)
    if args.verbose:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    if args.complete_check_cvmfs:
        is_build_complete = CVMFSPublishIsComplete(get_cvmfs_publish_records(args))
    else:
        is_build_complete = GalaxyHistoryIsBuildComplete(get_galaxy_history_names(args))

    split_options = SplitOptions()
    split_options.data_managers_path = args.data_managers_path
    split_options.merged_genomes_path = args.merged_genomes_path
    split_options.split_genomes_path = args.split_genomes_path
    split_options.is_build_complete = is_build_complete
    split_options.tool_id_mode = args.tool_id_mode

    filters = Filters()
    filters.build_id = args.filter_build_id
    filters.data_manager = args.filter_data_manager
    filters.stage = args.filter_stage
    split_options.filters = filters

    split_genomes(split_options)


if __name__ == "__main__":
    main()
