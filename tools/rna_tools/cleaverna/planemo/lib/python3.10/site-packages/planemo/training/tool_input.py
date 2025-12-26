"""Module contains code for the ToolInput class, dealing with the description of tool in workflow and XML."""

from planemo import templates
from planemo.io import info

INPUT_PARAM = """
>{{space}}- *"{{param_label}}"*: `{{param_value}}`
"""

INPUT_FILE_TEMPLATE = """
>{{space}}- {{ '{%' }} icon {{icon}} {{ '%}' }} *"{{input_name}}"*: {{input_value}}
"""

INPUT_SECTION = """
>{{space}}- In *"{{section_label}}"*:
"""

INPUT_ADD_REPEAT = """
>{{space}}- {{ '{%' }} icon param-repeat {{ '%}' }} *"Insert {{repeat_label}}"*
"""

SPACE = "    "


class ToolInput:
    """Class to describe a tool input / parameter and its value from a workflow."""

    def __init__(self, tool_inp_desc, wf_param_values, wf_steps, level, should_be_there=False, force_default=False):
        """Init an instance of ToolInput."""
        self.name = tool_inp_desc["name"]
        if "type" not in tool_inp_desc:
            raise ValueError(f"No type for the parameter {tool_inp_desc['name']}")
        self.type = tool_inp_desc["type"]
        self.tool_inp_desc = tool_inp_desc
        self.level = level
        self.wf_param_values = wf_param_values
        self.wf_steps = wf_steps
        self.formatted_desc = ""
        self.force_default = force_default

        if self.name not in self.wf_param_values:
            if not should_be_there:
                info(f"{self.name} not in workflow")
            else:
                raise ValueError(f"{self.name} not in workflow")
        else:
            self.wf_param_values = self.wf_param_values[self.name]

    def get_formatted_inputs(self):
        """Format the inputs of a step."""
        inputlist = ""
        inps = []
        if isinstance(self.wf_param_values, list):
            # multiple input (not collection)
            icon = "param-files"
            for i in self.wf_param_values:
                inps.append(f"`{i['output_name']}` {get_input_tool_name(i['id'], self.wf_steps)}")
        else:
            inp = self.wf_param_values
            if "id" in inp:
                # sinle input or collection
                inp_type = self.wf_steps[str(inp["id"])]["type"]
                if "collection" in inp_type:
                    icon = "param-collection"
                else:
                    icon = "param-file"
                inps = [f"`{inp['output_name']}` {get_input_tool_name(inp['id'], self.wf_steps)}"]
        if len(inps) > 0:
            inputlist += templates.render(
                INPUT_FILE_TEMPLATE,
                **{
                    "icon": icon,
                    "input_name": self.tool_inp_desc["label"],
                    "input_value": ", ".join(inps),
                    "space": SPACE * self.level,
                },
            )
        return inputlist

    def get_lower_param_desc(self):
        """Get the formatted description of the paramaters in the 'inputs' of the tool description."""
        sub_param_desc = ""
        for inp in self.tool_inp_desc["inputs"]:
            tool_inp = ToolInput(inp, self.wf_param_values, self.wf_steps, self.level + 1)
            sub_param_desc += tool_inp.get_formatted_desc()
        return sub_param_desc

    def get_formatted_section_desc(self):
        """Format the description (label and value) for parameters in a section."""
        section_paramlist = ""
        sub_param_desc = self.get_lower_param_desc()
        if sub_param_desc != "":
            section_paramlist += templates.render(
                INPUT_SECTION, **{"space": SPACE * self.level, "section_label": self.tool_inp_desc["title"]}
            )
            section_paramlist += sub_param_desc
        return section_paramlist

    def get_formatted_conditional_desc(self):
        """Format the description (label and value) for parameters in a conditional."""
        conditional_paramlist = ""
        # Get conditional parameter
        inp = ToolInput(
            self.tool_inp_desc["test_param"],
            self.wf_param_values,
            self.wf_steps,
            self.level,
            should_be_there=True,
            force_default=True,
        )
        conditional_paramlist = inp.get_formatted_desc()
        cond_param = inp.wf_param_values

        # Get parameters in the when and their values
        tmp_tool_inp_desc = self.tool_inp_desc
        for case in tmp_tool_inp_desc["cases"]:
            if case["value"] == cond_param and len(case["inputs"]) > 0:
                self.tool_inp_desc = case
                conditional_paramlist += self.get_lower_param_desc()
        self.tool_inp_desc = tmp_tool_inp_desc
        return conditional_paramlist

    def get_formatted_repeat_desc(self):
        """Format the description (label and value) for parameters in a repeat."""
        repeat_paramlist = ""
        if self.wf_param_values != "[]":
            tool_inp = {}
            for inp in self.tool_inp_desc["inputs"]:
                tool_inp.setdefault(inp["name"], inp)
            tmp_wf_param_values = self.wf_param_values
            cur_level = self.level
            for param in tmp_wf_param_values:
                self.wf_param_values = param
                self.level = cur_level + 1
                paramlist_in_repeat = self.get_lower_param_desc()
                if paramlist_in_repeat != "":
                    # add first click
                    repeat_paramlist += templates.render(
                        INPUT_ADD_REPEAT, **{"space": SPACE * (self.level), "repeat_label": self.tool_inp_desc["title"]}
                    )
                    repeat_paramlist += paramlist_in_repeat
                self.level = cur_level
            self.wf_param_values = tmp_wf_param_values

        repeat_desc = ""
        if repeat_paramlist != "":
            repeat_desc += (
                templates.render(
                    INPUT_SECTION, **{"space": SPACE * self.level, "section_label": self.tool_inp_desc["title"]}
                )
                + repeat_paramlist
            )
        return repeat_desc

    def get_formatted_other_param_desc(self):
        """Get value of a 'simple' parameter if different from the default value, None otherwise."""
        param_value = None
        if self.tool_inp_desc["value"] == self.wf_param_values and not self.force_default:
            param_value = None
        elif self.type == "boolean":
            if bool(self.tool_inp_desc["value"]) == self.wf_param_values:
                param_value = None
            else:
                param_value = "Yes" if self.wf_param_values else "No"
        elif self.type == "select":
            param_values = []
            for opt in self.tool_inp_desc["options"]:
                if opt[1] == self.wf_param_values:
                    param_values.append(opt[0])
            param_value = ", ".join(param_values)
        elif self.type == "data_column":
            param_value = f"c{self.wf_param_values}"
        else:
            param_value = self.wf_param_values

        param_desc = ""
        if param_value is not None:
            param_desc = templates.render(
                INPUT_PARAM,
                **{"space": SPACE * self.level, "param_label": self.tool_inp_desc["label"], "param_value": param_value},
            )
        return param_desc

    def get_formatted_desc(self):
        """Get the formatted description (ready for hands-on tutorial) of the parameter."""
        if self.wf_param_values:
            if self.type == "data" or self.type == "data_collection":
                self.formatted_desc += self.get_formatted_inputs()
            elif self.type == "section":
                self.formatted_desc += self.get_formatted_section_desc()
            elif self.type == "conditional":
                self.formatted_desc += self.get_formatted_conditional_desc()
            elif self.type == "repeat":
                self.formatted_desc += self.get_formatted_repeat_desc()
            else:
                self.formatted_desc += self.get_formatted_other_param_desc()
        return self.formatted_desc


def get_input_tool_name(step_id, steps):
    """Get the string with the name of the tool that generated an input."""
    inp_provenance = ""
    inp_prov_id = str(step_id)
    if inp_prov_id in steps:
        name = steps[inp_prov_id]["name"]
        if "Input dataset" in name:
            inp_provenance = f"({name})"
        else:
            inp_provenance = f"(output of **{name}** {{% icon tool %}})"
    return inp_provenance


def get_empty_input():
    """Get the string for an empty input."""
    return templates.render(
        INPUT_FILE_TEMPLATE,
        **{"space": 1 * SPACE, "icon": "param-file", "input_name": "Input file", "input_value": "File"},
    )


def get_empty_param():
    """Get the string for an empty param."""
    return templates.render(INPUT_PARAM, **{"space": 1 * SPACE, "param_label": "Parameter", "param_value": "a value"})
