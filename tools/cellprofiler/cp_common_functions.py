INDENTATION = "    "
LINE_NUM_MODULES = 4

def get_json_value(input_params, keys_path):
    """Returns the value in keys_path or empty string if the field does not exist"""
    if not isinstance(input_params, dict):
        return ""
    key_list = keys_path.split(".")
    try:
        value = input_params[key_list[0]]
        for k in key_list[1:]:
            value = value[k]
        return value
    except KeyError:
        return ""


def concat_conditional(a, b):
    if a == "" or b == "":
        return ""
    else:
        return f"{a}_{b}"

        
def get_total_number_of_modules(pipeline_lines):
    """Gets the number of modules from the header of the previous pipeline"""
    number_of_modules = pipeline_lines[LINE_NUM_MODULES].strip().split(':')[1]
    return int(number_of_modules)


def get_pipeline_lines(input_pipeline):
    """Returns a list with the lines in the .cppipe file"""
    with open(input_pipeline) as f:
        lines = f.readlines()
    return lines


def update_module_count(pipeline_lines, count):
    """Updates the number of modules in the .cppipe header"""
    module_count_entry = pipeline_lines[LINE_NUM_MODULES].strip().split(':')[0]
    pipeline_lines[4] = f"{module_count_entry}:{count}\n"
    return pipeline_lines


def write_pipeline(filename, lines_pipeline):
    """Writes the lines composing the pipeline into a file"""
    with open(filename, "w") as f:
        f.writelines(lines_pipeline)


def build_header(module_name, module_number):
    """Creates the first line of a module given the name and module number"""
    result = "|".join([f"{module_name}:[module_num:{module_number}",
                       "svn_version:\\'Unknown\\'",
                       "variable_revision_number:4",
                       "show_window:False",
                       "notes:\\x5B\\x5D",
                       "batch_state:array(\\x5B\\x5D, dtype=uint8)",
                       "enabled:True",
                       "wants_pause:False]\n"])
    return result
