#!/usr/bin/env python

import argparse
import json

from cp_common_functions import get_json_value
from cp_common_functions import get_pipeline_lines
from cp_common_functions import get_total_number_of_modules
from cp_common_functions import INDENTATION
from cp_common_functions import update_module_count
from cp_common_functions import write_pipeline

MODULE_NAME = "ColorToGray"
OUTPUT_FILENAME = "output.cppipe"


def build_ctg_header(module_name, module_number):
    """Creates the first line of a module given the name and module number"""
    result = "|".join(
        [
            f"{module_name}:[module_num:{module_number}",
            "svn_version:\\'Unknown\\'",
            "variable_revision_number:4",
            "show_window:True",
            "notes:\\x5B\\'Convert the color image to grayscale.\\'\\x5D",
            "batch_state:array(\\x5B\\x5D, dtype=uint8)",
            "enabled:True",
            "wants_pause:False]\n",
        ]
    )
    return result


def build_main_block(input_params):
    """Creates the main block of the CP pipeline with the parameters that don't depend on conditional choices"""
    result = INDENTATION.join(
        [
            f"{INDENTATION}Select the input image:{get_json_value(input_params,'name_input_image')}\n",
            f"Conversion method:{get_json_value(input_params,'con_conversion_method.conversion_method')}\n",
            f"Image type:{get_json_value(input_params,'con_conversion_method.con_image_type.image_type')}\n",
        ]
    )

    conversion_method = get_json_value(
        input_params, "con_conversion_method.conversion_method"
    )

    image_type = get_json_value(
        input_params, "con_conversion_method.con_image_type.image_type"
    )
    rgb_comb_name_out = "OrigGray"
    combine_w_red = 1.0
    combine_w_green = 1.0
    combine_w_blue = 1.0

    split_red = "Yes"
    split_green = "Yes"
    split_blue = "Yes"
    red_output_name = "OrigRed"
    green_output_name = "OrigGreen"
    blue_output_name = "OrigBlue"

    split_hue = "Yes"
    split_saturation = "Yes"
    split_value = "Yes"
    hue_output_name = "OrigHue"
    saturation_output_name = "OrigSaturation"
    value_output_name = "OrigValue"

    channel_count = 1
    if conversion_method == "Combine":
        if image_type == "RGB" or image_type == "HSV":
            rgb_comb_name_out = get_json_value(
                input_params, "con_conversion_method.name_output_image"
            )
            combine_w_red = get_json_value(
                input_params, "con_conversion_method.con_image_type.weight_red_channel"
            )
            combine_w_green = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.weight_green_channel",
            )
            combine_w_blue = get_json_value(
                input_params, "con_conversion_method.con_image_type.weight_blue_channel"
            )
    elif conversion_method == "Split":
        if image_type == "RGB":
            split_red = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_red.yes_no",
            )
            red_output_name = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_red.name_output_image",
            )
            split_green = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_green.yes_no",
            )
            green_output_name = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_green.name_output_image",
            )
            split_blue = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_blue.yes_no",
            )
            blue_output_name = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_blue.name_output_image",
            )
        elif image_type == "HSV":
            split_hue = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_hue.yes_no",
            )
            hue_output_name = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_hue.name_output_image",
            )
            split_saturation = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_saturation.yes_no",
            )
            saturation_output_name = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_saturation.name_output_image",
            )
            split_value = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_value.yes_no",
            )
            value_output_name = get_json_value(
                input_params,
                "con_conversion_method.con_image_type.con_convert_value.name_output_image",
            )

    result += INDENTATION.join(
        [
            f"{INDENTATION}Name the output image:{rgb_comb_name_out}\n",
            f"Relative weight of the red channel:{str(combine_w_red)}\n",
            f"Relative weight of the green channel:{str(combine_w_green)}\n",
            f"Relative weight of the blue channel:{str(combine_w_blue)}\n",
            f"Convert red to gray?:{split_red}\n",
            f"Name the output image:{red_output_name}\n",
            f"Convert green to gray?:{split_green}\n",
            f"Name the output image:{green_output_name}\n",
            f"Convert blue to gray?:{split_blue}\n",
            f"Name the output image:{blue_output_name}\n",
            f"Convert hue to gray?:{split_hue}\n",
            f"Name the output image:{hue_output_name}\n",
            f"Convert saturation to gray?:{split_saturation}\n",
            f"Name the output image:{saturation_output_name}\n",
            f"Convert value to gray?:{split_value}\n",
            f"Name the output image:{value_output_name}\n",
        ]
    )

    channel_count = 1
    if image_type == "Channels":
        channels = input_params["con_conversion_method"]["con_image_type"][
            "rpt_channel"
        ]
        channel_count = len(channels)
        result += INDENTATION.join([f"{INDENTATION}Channel count:{channel_count}\n"])

        for ch in channels:
            rel_weight_ch = 1.0
            image_name = "Channel1"
            if conversion_method == "Combine":
                rel_weight_ch = get_json_value(ch, "weight_of_channel")
            else:
                image_name = get_json_value(ch, "image_name")
            result += INDENTATION.join(
                [
                    f"{INDENTATION}Channel number:{get_json_value(ch,'channel_no')}\n",
                    f"Relative weight of the channel:{str(rel_weight_ch)}\n",
                    f"Image name:{image_name}\n",
                ]
            )
    else:
        result += INDENTATION.join(
            [
                f"{INDENTATION}Channel count:{channel_count}\n",
                "Channel number:Red\\x3A 1\n",
                "Relative weight of the channel:1.0\n",
                "Image name:Channel1\n",
            ]
        )
    result = result.rstrip("\n")
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
