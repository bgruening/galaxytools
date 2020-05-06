import json
import sys
import os

FOURSPACES = "    "

input_json_path = sys.argv[1]

params = json.load(open(input_json_path, "r"))


def write_images():
    filter_images = params['images']['filter_images']

    _str = "\nImages:[module_num:1|svn_version:\\'Unknown\\'|variable_revision_number:2|show_window:False|notes:\\x5B\\'To begin creating your project, use the Images module to compile a list of files and/or folders that you want to analyze. You can also specify a set of rules to include only the desired files in your selected folders.\\'\x5D|batch_state:array(\x5B\x5D, dtype=uint8)|enabled:True|wants_pause:False]\n"
    _str += FOURSPACES+":\n"
    _str += FOURSPACES + "Filter images?:%s\n" % filter_images
    _str += FOURSPACES + "Select the rule criteria:and (extension does isimage) (directory doesnot startwith \".\")\n"

    return _str


def write_metadata():
    metadata_extraction = params['metadata']['con_metadata_extraction']
    extract = metadata_extraction['extract']

    if 'extraction_method' in metadata_extraction:
        method_count = len(metadata_extraction['extraction_method'])
    else:
        method_count = 1

    _str = "\nMetadata:[module_num:2|svn_version:\\'Unknown\\'|variable_revision_number:4|show_window:False|notes:\\x5B\\'The Metadata module optionally allows you to extract information describing your images (i.e, metadata) which will be stored along with your measurements. This information can be contained in the file name and/or location, or in an external file.\\'\x5D|batch_state:array(\x5B\x5D, dtype=uint8)|enabled:True|wants_pause:False]\n"
    _str += FOURSPACES + "Extract metadata?:%s\n" % extract

    if extract == "No":
        _str += FOURSPACES + "Metadata data type:Text\n"
        _str += FOURSPACES + "Metadata types:{}\n"
        _str += FOURSPACES + "Extraction method count:%d\n" % method_count
        _str += FOURSPACES + "Metadata extraction method:Extract from file/folder names\n"
        _str += FOURSPACES + "Regular expression to extract from file name:^(?P<Plate>.*)_(?P<Well>\x5BA-P\x5D\x5B0-9\x5D{2})_s(?P<Site>\x5B0-9\x5D)_w(?P<ChannelNumber>\x5B0-9\x5D)\n"
        _str += FOURSPACES + "Regular expression to extract from folder name:(?P<Date>\x5B0-9\x5D{4}_\x5B0-9\x5D{2}_\x5B0-9\x5D{2})$\n"
        _str += FOURSPACES + "Extract metadata from:All images\n"
        _str += FOURSPACES + "Select the filtering criteria:and (file does contain \"\")\n"
        _str += FOURSPACES + "Metadata file location:\n"
        _str += FOURSPACES + "Match file and image metadata:\x5B\x5D\n"
        _str += FOURSPACES + "Use case insensitive matching?:No\n"
    else:
        _str += FOURSPACES + "Metadata data type:Text\n"  #default Text,not possible to select in Galaxy
        _str += FOURSPACES + "Metadata types:{}\n"
        _str += FOURSPACES + "Extraction method count:%d\n" % method_count

        for methods in metadata_extraction["extraction_method"]:
            _str += FOURSPACES + "Metadata extraction method:%s\n" % methods["metadata_extraction_method"]
            _str += FOURSPACES + "Metadata source:%s\n" % methods["con_metadata_source"]["metadata_source"]

            if "file_name_regex" in methods["con_metadata_source"]:
                file_regex = methods["con_metadata_source"]["file_name_regex"]
                folder_regex = "(?P<folderField1>.*)"
            elif "folder_name_regex" in methods["con_metadata_source"]:
                file_regex = "(?P<field1>.*)_(?P<field2>[a-zA-Z0-9]+)"
                folder_regex = methods["con_metadata_source"]["folder_name_regex"]
            else:
                # default regex
                file_regex = "(?P<field1>.*)_(?P<field2>[a-zA-Z0-9]+)"
                folder_regex = "(?P<folderField1>.*)"

            _str += FOURSPACES + "Regular expression to extract from file name:%s\n" % file_regex
            _str += FOURSPACES + "Regular expression to extract from folder name:%s\n" % folder_regex

            _str += FOURSPACES + "Extract metadata from:%s\n" % methods["con_metadata_extract_from"]["extract_metadata_from"]

            if methods["con_metadata_extract_from"]["extract_metadata_from"] == "Images matching a rule":
                rule_str =""
                for r in methods["con_metadata_extract_from"]["r_match"]:
                    if r["con_match"]["rule_type"] == "extension":
                        rule_str += " (" + r["con_match"]["rule_type"] + " " + r["con_match"]["operator"] + " " + \
                                    r["con_match"]["match_type"]+")"
                    else:
                        rule_str +=" (" + r["con_match"]["rule_type"] + " " + r["con_match"]["operator"] + " " +\
                                   r["con_match"]["contain"] + " \"" + r["con_match"]["match_value"] +"\")"


                _str += FOURSPACES + "Select the filtering criteria:" + methods["con_metadata_extract_from"]["match_all_any"] + rule_str +"\n"
            else:
                _str += FOURSPACES + "Select the filtering criteria:and (file does contain \"\")\n" #this line is required even if it's not used

            _str += FOURSPACES + "Metadata file location:\n"
            _str += FOURSPACES + "Match file and image metadata:\x5B\x5D\n"
            _str += FOURSPACES + "Use case insensitive matching?:No\n"

    return _str


def write_nameandtypes():
    nameandtypes = params['nameandtypes']
    assign_a_name = nameandtypes['con_assign_a_name_to']['assign_a_name_to']

    if "con_select_image_type" in nameandtypes['con_assign_a_name_to']:
        con_set_intensity = nameandtypes['con_assign_a_name_to']['con_select_image_type']['con_set_intensity']
        max_intensity = con_set_intensity['maximum_intensity'] if "maximum_intensity" in con_set_intensity else 255.0
    else:
        max_intensity = 255.0

    pixel_space = nameandtypes['pixel_space']

    rule_count = len(nameandtypes['con_assign_a_name_to']['r_match_rule']) if "r_match_rule" in nameandtypes['con_assign_a_name_to'] else 1

    process_3d = nameandtypes['pixel_space']['process_3d']
    x_spacing = 1.0 if "x_spacing" not in pixel_space else pixel_space["x_spacing"]
    y_spacing = 1.0 if "y_spacing" not in pixel_space else pixel_space["y_spacing"]
    z_spacing = 1.0 if "z_spacing" not in pixel_space else pixel_space["z_spacing"]

    _str = "\nNamesAndTypes:[module_num:3|svn_version:\\'Unknown\\'|variable_revision_number:8|show_window:False|notes:\\x5B\\'The NamesAndTypes module allows you to assign a meaningful name to each image by which other modules will refer to it.\\'\\x5D|batch_state:array(\\x5B\\x5D, dtype=uint8)|enabled:True|wants_pause:False]\n"

    _str += FOURSPACES + "Assign a name to:%s\n" % assign_a_name

    if assign_a_name == "All images":
        _str += FOURSPACES + "Select the image type:%s\n" % nameandtypes['con_assign_a_name_to']['con_select_image_type']['select_image_type']
        _str += FOURSPACES + "Name to assign these images:%s\n" % nameandtypes['con_assign_a_name_to']['name_to_assign']
        _str += FOURSPACES + "Match metadata:[]\n"

        _str += FOURSPACES + "Image set matching method:Order\n"
        _str += FOURSPACES + "Set intensity range from:%s\n" % con_set_intensity['set_intensity_range_from']
        _str += FOURSPACES + "Assignments count:%s\n" % rule_count
        _str += FOURSPACES + "Single images count:0\n"
        _str += FOURSPACES + "Maximum intensity:%.1f\n" % max_intensity
        _str += FOURSPACES + "Process as 3D?:%s\n" % process_3d

    else:
        #the below lines are not relevant to "images matching rules", but needed in pipeline file
        _str += FOURSPACES + "Select the image type:Grayscale image\n"
        _str += FOURSPACES + "Name to assign these images:DNA\n"
        _str += FOURSPACES + "Match metadata:[]\n"

        _str += FOURSPACES + "Image set matching method:%s\n" % nameandtypes['con_assign_a_name_to']['matching_method']
        _str += FOURSPACES + "Set intensity range from:Image metadata\n"
        _str += FOURSPACES + "Assignments count:%d\n" % rule_count
        _str += FOURSPACES + "Single images count:0\n"
        _str += FOURSPACES + "Maximum intensity:%.1f\n" % max_intensity
        _str += FOURSPACES + "Process as 3D?:%s\n" % process_3d

    _str += FOURSPACES + "Relative pixel spacing in X:%.1f\n" % x_spacing
    _str += FOURSPACES + "Relative pixel spacing in Y:%.1f\n" % y_spacing
    _str += FOURSPACES + "Relative pixel spacing in Z:%.1f\n" % z_spacing

    if assign_a_name == "Images matching rules":
        for rule in nameandtypes["con_assign_a_name_to"]["r_match_rule"]:

            rule_str = ""
            if len(rule["r_match"]) >0 :
                for r in rule["r_match"]:
                        if r["con_match"]["rule_type"] == "file" or r["con_match"]["rule_type"] == "directory":
                            rule_str += " (" + r["con_match"]["rule_type"] + " "+r["con_match"]["operator"]+" "+\
                                        r["con_match"]["contain"]+" \"" + r["con_match"]["match_value"] +"\")"
                        else:
                            rule_str += " ("+ r["con_match"]["rule_type"] + " " + r["con_match"]["operator"] + " " + \
                                        r["con_match"]["match_type"] + ")"
            else:
                rule_str = " (file does contain \"\")"  #need to have a value even if it is not used

            _str += FOURSPACES + "Select the rule criteria:" + rule["match_all_any"] + rule_str +"\n"

            img_or_obj = rule["con_select_image_type"]["select_image_type"]

            if img_or_obj == "Objects":
                _str += FOURSPACES + "Name to assign these images:DNA\n"
                _str += FOURSPACES + "Name to assign these objects:%s\n" % rule["con_select_image_type"]["name_to_assign"]
            else:
                _str += FOURSPACES + "Name to assign these images:%s\n" % rule["con_select_image_type"]["name_to_assign"]
                _str += FOURSPACES + "Name to assign these objects:Cell\n"

            _str += FOURSPACES + "Select the image type:%s\n" % img_or_obj


            intensity_range="Image metadata" #default value
            if img_or_obj == "Grayscale image" or img_or_obj == "Color image":
                intensity_range = rule["con_select_image_type"]["con_set_intensity"]["set_intensity_range_from"]

            _str += FOURSPACES + "Set intensity range from:%s\n" % intensity_range

            if intensity_range == "Manual":
                _str += FOURSPACES + "Maximum intensity:%s\n" % rule["con_select_image_type"]["con_set_intensity"]["maximum_intensity"]
            else:
                _str += FOURSPACES + "Maximum intensity:255.0\n"


    return _str


def write_groups():
    groups = params['groups']

    _str = "\nGroups:[module_num:4|svn_version:\\'Unknown\\'|variable_revision_number:2|show_window:False|notes:\\x5B\\\'The Groups module optionally allows you to split your list of images into image subsets (groups) which will be processed independently of each other. Examples of groupings include screening batches, microtiter plates, time-lapse movies, etc.\\'\\x5D|batch_state:array(\\x5B\\x5D, dtype=uint8)|enabled:True|wants_pause:False]\n"

    group_images =  groups["con_groups"]["group_images"]

    _str += FOURSPACES + "Do you want to group your images?:%s\n" % group_images
    _str += FOURSPACES + "grouping metadata count:1\n"

    if group_images == "Yes":
        _str += FOURSPACES + "Metadata category:%s\n" % groups["con_groups"]["group_category"]
    else:
        _str += FOURSPACES + "Metadata category:None\n"

    return _str


with open("output.cppipe", "w") as f:
    headers = ["CellProfiler Pipeline: http://www.cellprofiler.org\n",
               "Version:3\n",
               "DateRevision:319\n",
               "GitHash:\n",
               "ModuleCount:4\n",
               "HasImagePlaneDetails:False",
               "\n"]

    f.writelines(headers)

    img_str = write_images()
    metadata_str = write_metadata()
    nameandtypes_str = write_nameandtypes()
    groups_str = write_groups()

    output_str = img_str + metadata_str + nameandtypes_str + groups_str

    f.write(output_str)
    f.close()