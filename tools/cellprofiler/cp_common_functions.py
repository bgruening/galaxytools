INDENTATION = "    "

def get_field(dictionary, path, default=None):
    keys = path.split("/")
    try:
        temp = dictionary[keys[0]]
        for k in keys[1:]:
            temp = dictionary[k]
        return(temp)
    except KeyError: pass


def get_total_number_of_modules(pipeline_lines):
    number_of_modules = pipeline_lines[4].strip().split(':')[1]
    return(int(number_of_modules))


def get_pipeline_lines(input_pipeline):
    with open(input_pipeline) as f:
        lines = f.readlines()
    return(lines)


def update_module_count(pipeline_lines, current_module_num):
    module_count_entry = pipeline_lines[4].strip().split(':')[0]
    pipeline_lines[4] = f"{module_count_entry}:{current_module_num}\n"
    return(pipeline_lines)
        

def write_pipeline(filename, lines_pipeline):
    with open(filename, "w") as f:
        f.writelines(lines_pipeline)


def build_header(module_number, module_name):
    result = "|".join([f"{module_name}:[module_num:{module_number}",
                       "svn_version:\\'Unknown\\'",
                       "variable_revision_number:4",
                       "show_window:False",
                       "notes:\\x5B\\x5D",
                       "batch_state:array(\\x5B\\x5D, dtype=uint8)",
                       "enabled:True",
                       "wants_pause:False]\n"])
    return(result)
