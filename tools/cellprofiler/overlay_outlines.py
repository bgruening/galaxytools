#!/usr/bin/env python

import argparse
import json

from cp_common_functions import (
    INDENTATION,
    get_json_value,
    get_pipeline_lines,
    get_total_number_of_modules,
    update_module_count,
    write_pipeline,
)

MODULE_NAME = "OverlayOutlines"
OUTPUT_FILENAME = "output.cppipe"


def build_ctg_header(module_name, module_number):
    """Creates the first line of a module given the name and module number"""
    result = "|".join(
        [
            f"{module_name}:[module_num:{module_number}",
            "svn_version:\\'Unknown\\'",
            "variable_revision_number:4",
            "show_window:True",
            "notes:\\x5B\\'Overlay the embryo outlines on the grayscale image.\\'\\x5D",
            "batch_state:array(\\x5B\\x5D, dtype=uint8)",
            "enabled:True",
            "wants_pause:False]\n",
        ]
    )
    return result


def build_main_block(input_params):
    result = f"{INDENTATION}Display outlines on a blank image?:{get_json_value(input_params,'con_blank_img.blank_img')}\n"

    on_blank = get_json_value(input_params, "con_blank_img.blank_img")
    # defaults
    img_to_display = "None"
    display_mode = get_json_value(
        input_params, "con_blank_img.con_display_mode.display_mode"
    )
    method_brightness = "Max of image"
    howto = get_json_value(input_params, "howto_outline")
    outline_color = "#FF0000"
    obj_to_display = "None"

    name_output_img = get_json_value(input_params, "name_output_image")

    if on_blank == "No":
        img_to_display = get_json_value(input_params, "con_blank_img.image_to_outline")

    result += INDENTATION.join(
        [
            f"{INDENTATION}Select image on which to display outlines:{img_to_display}\n",
            f"Name the output image:{name_output_img}\n",
            f"Outline display mode:{display_mode}\n",
        ]
    )

    if on_blank == "No" and display_mode == "Grayscale":
        method_brightness = get_json_value(
            input_params, "con_blank_img.con_display_mode.method_brightness"
        )

    result += INDENTATION.join(
        [
            f"{INDENTATION}Select method to determine brightness of outlines:{method_brightness}\n",
            f"How to outline:{howto}\n",
        ]
    )

    obj_outline_str = ""

    if display_mode == "Color":
        for obj in input_params["con_blank_img"]["con_display_mode"][
            "rpt_obj_to_display"
        ]:
            outline_color = get_json_value(obj, "outline_color")
            obj_to_display = get_json_value(obj, "obj_to_display")
            obj_outline_str += INDENTATION.join(
                [
                    f"{INDENTATION}Select outline color:{outline_color}\n",
                    f"Select objects to display:{obj_to_display}\n",
                ]
            )
    else:  # grayscale
        for obj in input_params["con_blank_img"]["con_display_mode"][
            "rpt_obj_to_display"
        ]:
            obj_to_display = get_json_value(obj, "obj_to_display")
            obj_outline_str += INDENTATION.join(
                [
                    f"{INDENTATION}Select outline color:{outline_color}\n",
                    f"Select objects to display:{obj_to_display}\n",
                ]
            )
    obj_outline_str = obj_outline_str.rstrip("\n")
    result += f"{obj_outline_str}"

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pipeline", help="CellProfiler pipeline")
    parser.add_argument("-i", "--inputs", help="JSON inputs from Galaxy")
    args = parser.parse_args()

    pipeline_lines = get_pipeline_lines(args.pipeline)
    inputs_galaxy = json.load(open(args.inputs, "r"))

    current_module_num = get_total_number_of_modules(pipeline_lines)
    current_module_num += 1
    pipeline_lines = update_module_count(pipeline_lines, current_module_num)

    header_block = build_ctg_header(MODULE_NAME, current_module_num)
    main_block = build_main_block(inputs_galaxy)

    module_pipeline = f"\n{header_block}{main_block}\n"
    pipeline_lines.append(module_pipeline)

    write_pipeline(OUTPUT_FILENAME, pipeline_lines)
