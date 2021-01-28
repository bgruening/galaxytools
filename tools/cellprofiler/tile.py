#!/usr/bin/env python

import json
import argparse
from cp_common_functions import * # noqa

MODULE_NAME = "Tile"
OUTPUT_FILENAME = "output.cppipe"


def build_header(module_name, module_number):
    result = "|".join([f"{module_name}:[module_num:{module_number}",
                       "svn_version:\\'Unknown\\'",
                       "variable_revision_number:1",
                       "show_window:True",
                       "notes:\\x5B\\'Tile the original color image, the outlined image and the image of tracked labels together.\\'\\x5D",
                       "batch_state:array(\\x5B\\x5D, dtype=uint8)",
                       "enabled:True",
                       "wants_pause:False]\n"])
    return result


def build_main_block(input_params):
    result = INDENTATION.join([f"{INDENTATION}Select an input image:{get_json_value(input_params,'input_image')}\n",
                               f"Name the output image:{get_json_value(input_params,'output_image_name')}\n",
                               f"Tile assembly method:{get_json_value(input_params,'con_assembly_method.assembly_method')}\n"
                               ])

    calc_rows = get_json_value(input_params, 'con_assembly_method.con_calc_no_row.calc_no_row')
    no_of_rows = 8

    calc_cols = get_json_value(input_params, 'con_assembly_method.con_calc_no_cols.calc_no_cols')
    no_of_cols = 12

    if calc_rows == "No":
        no_of_rows = get_json_value(input_params, 'con_assembly_method.con_calc_no_row.no_of_row')

    if calc_cols == "No":
        no_of_cols = get_json_value(input_params, 'con_assembly_method.con_calc_no_cols.no_of_cols')

    corner_to_begin = get_json_value(input_params, 'con_assembly_method.corner_to_begin')
    direction_tiling = get_json_value(input_params, 'con_assembly_method.direction')
    meander = get_json_value(input_params, 'con_assembly_method.meander_mode')

    assembly_method = get_json_value(input_params, 'con_assembly_method.assembly_method')

    result += INDENTATION.join(
        [f"{INDENTATION}Final number of rows:{str(no_of_rows)}\n",
         f"Final number of columns:{str(no_of_cols)}\n",
         f"Image corner to begin tiling:{corner_to_begin}\n",
         f"Direction to begin tiling:{direction_tiling}\n",
         f"Use meander mode?:{meander}\n",
         f"Automatically calculate number of rows?:{calc_rows}\n",
         f"Automatically calculate number of columns?:{calc_cols}\n"
         ])

    if assembly_method == "Within cycles":
        additionals = input_params['con_assembly_method']['rpt_additional_image']

        for img in additionals:
            result += INDENTATION.join(
                [f"{INDENTATION}Select an additional image to tile:{get_json_value(img, 'additional_img')}\n"
                 ])

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-p', '--pipeline',
        help='CellProfiler pipeline'
    )
    parser.add_argument(
        '-i', '--inputs',
        help='JSON inputs from Galaxy'
    )
    args = parser.parse_args()

    pipeline_lines = get_pipeline_lines(args.pipeline)
    inputs_galaxy = json.load(open(args.inputs, "r"))

    current_module_num = get_total_number_of_modules(pipeline_lines)
    current_module_num += 1
    pipeline_lines = update_module_count(pipeline_lines, current_module_num)

    header_block = build_header(MODULE_NAME, current_module_num)
    main_block = build_main_block(inputs_galaxy)

    module_pipeline = f"\n{header_block}{main_block}\n"
    pipeline_lines.append(module_pipeline)

    write_pipeline(OUTPUT_FILENAME, pipeline_lines)
