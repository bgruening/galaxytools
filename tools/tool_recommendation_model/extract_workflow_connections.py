"""
Extract workflow paths from the tabular file containing
input and output tools
"""

import csv
import random

import utils


class ExtractWorkflowConnections:
    def __init__(self):
        """ Init method. """

    def collect_standard_connections(self, row):
        published = row[8]
        deleted = row[9]
        has_errors = row[10]
        if published == "t" and deleted == "f" and has_errors == "f":
            return True
        return False

    def read_tabular_file(self, raw_file_path):
        """
        Read tabular file and extract workflow connections
        """
        print("Reading workflows...")
        workflows = {}
        workflow_paths_dup = ""
        workflow_parents = dict()
        workflow_paths = list()
        unique_paths = dict()
        standard_connections = dict()
        with open(raw_file_path, "rt") as workflow_connections_file:
            workflow_connections = csv.reader(workflow_connections_file, delimiter="\t")
            for index, row in enumerate(workflow_connections):
                wf_id = str(row[0])
                in_tool = row[3].strip()
                out_tool = row[6].strip()
                if wf_id not in workflows:
                    workflows[wf_id] = list()
                if out_tool and in_tool and out_tool != in_tool:
                    workflows[wf_id].append((out_tool, in_tool))
                    qc = self.collect_standard_connections(row)
                    if qc:
                        i_t = utils.format_tool_id(in_tool)
                        o_t = utils.format_tool_id(out_tool)
                        if i_t not in standard_connections:
                            standard_connections[i_t] = list()
                        if o_t not in standard_connections[i_t]:
                            standard_connections[i_t].append(o_t)
        print("Processing workflows...")
        wf_ctr = 0
        for wf_id in workflows:
            wf_ctr += 1
            workflow_parents[wf_id] = self.read_workflow(wf_id, workflows[wf_id])

        for wf_id in workflow_parents:
            flow_paths = list()
            parents_graph = workflow_parents[wf_id]
            roots, leaves = self.get_roots_leaves(parents_graph)
            for root in roots:
                for leaf in leaves:
                    paths = self.find_tool_paths_workflow(parents_graph, root, leaf)
                    # reverse the paths as they are computed from leaves to roots leaf
                    paths = [tool_path for tool_path in paths]
                    if len(paths) > 0:
                        flow_paths.extend(paths)
            workflow_paths.extend(flow_paths)
        print("Workflows processed: %d" % wf_ctr)

        # remove slashes from the tool ids
        wf_paths_no_slash = list()
        for path in workflow_paths:
            path_no_slash = [utils.format_tool_id(tool_id) for tool_id in path]
            wf_paths_no_slash.append(path_no_slash)

        # collect duplicate paths
        for path in wf_paths_no_slash:
            workflow_paths_dup += ",".join(path) + "\n"

        # collect unique paths
        unique_paths = list(workflow_paths_dup.split("\n"))
        unique_paths = list(filter(None, unique_paths))
        random.shuffle(unique_paths)
        no_dup_paths = list(set(unique_paths))

        print("Finding compatible next tools...")
        compatible_next_tools = self.set_compatible_next_tools(no_dup_paths)
        return unique_paths, compatible_next_tools, standard_connections

    def set_compatible_next_tools(self, workflow_paths):
        """
        Find next tools for each tool
        """
        next_tools = dict()
        for path in workflow_paths:
            path_split = path.split(",")
            for window in range(0, len(path_split) - 1):
                current_next_tools = path_split[window: window + 2]
                current_tool = current_next_tools[0]
                next_tool = current_next_tools[1]
                try:
                    next_tools[current_tool] += "," + next_tool
                except Exception:
                    next_tools[current_tool] = next_tool
        for tool in next_tools:
            next_tools[tool] = ",".join(list(set(next_tools[tool].split(","))))
        return next_tools

    def read_workflow(self, wf_id, workflow_rows):
        """
        Read all connections for a workflow
        """
        tool_parents = dict()
        for connection in workflow_rows:
            in_tool = connection[0]
            out_tool = connection[1]
            if out_tool not in tool_parents:
                tool_parents[out_tool] = list()
            if in_tool not in tool_parents[out_tool]:
                tool_parents[out_tool].append(in_tool)
        return tool_parents

    def get_roots_leaves(self, graph):
        roots = list()
        leaves = list()
        all_parents = list()
        for item in graph:
            all_parents.extend(graph[item])
        all_parents = list(set(all_parents))
        children = graph.keys()
        roots = list(set(all_parents).difference(set(children)))
        leaves = list(set(children).difference(set(all_parents)))
        return roots, leaves

    def find_tool_paths_workflow(self, graph, start, end, path=[]):
        path = path + [end]
        if start == end:
            return [path]
        path_list = list()
        if end in graph:
            for node in graph[end]:
                if node not in path:
                    new_tools_paths = self.find_tool_paths_workflow(
                        graph, start, node, path
                    )
                    for tool_path in new_tools_paths:
                        path_list.append(tool_path)
        return path_list
