import json
import sys
import os
import argparse
from cp_common_functions import *

MODULE_NAME = "ImageMath"
OUTPUT_FILENAME = "output.cppipe"

operator_map = {
    "add": "Add",
    "subtract": "Subtract",
    "multiply": "Multiply",
    "divide": "Divide",
    "average": "Average",
    "minimum": "Minimum",
    "maximum": "Maximum",
    "invert": "Invert",
    "log_2": "Log transform (base 2)",
    "log_legacy": "Log transform (legacy)",
    "and": "And",
    "or": "Or",
    "not": "Not",
    "equals": "Equals"
}


def build_main_block(input_params):
    operation = operator_map[input_params['operation']['operation']]
    result = INDENTATION.join(
        [f"{INDENTATION}Operation:{operation.encode('utf-16')}\n",
         f"Raise the power of the result by:{input_params['operation']['op_results']['raise_the_power_of_the_result_by']}\n",
         f"Multiply the result by:{input_params['operation']['op_results']['multiply_the_result_by']}\n",
         f"Add to result:{input_params['operation']['op_results']['add_to_result']}\n",
         f"Set values less than 0 equal to 0?:{input_params['operation']['op_results']['set_values_less_than_0_equal_to_0']}\n",
         f"Set values greater than 1 equal to 1?:{input_params['operation']['op_results']['set_values_greater_than_1_equal_to_1']}\n",
         f"Ignore the image masks?:{input_params['ignore_the_image_masks']}\n",
         f"Name the output image:{input_params['name_output_image']}"
         ])
    return(result)


def build_first_image_block(input_params):
    result = INDENTATION.join(
        [f"{INDENTATION}Image or measurement?:{input_params['image_or_measurement']['image_or_measurement']}\n",
         f"Select the first image:{input_params['image_or_measurement']['select_the_first_image']}\n",
         f"Multiply the first image by:{input_params['multiply_the_first_image_by']}\n",
         f"Measurement:{get_field(input_params,'image_or_measurement/measurement')}"
         ])
    return(result)


def build_second_image_block(input_params):
    result = INDENTATION.join(
        [f"{INDENTATION}Image or measurement?:{input_params['image_or_measurement']['image_or_measurement']}\n",
         f"Select the second image:{input_params['image_or_measurement']['select_the_second_image']}\n",
         f"Multiply the second image by:{get_field(input_params,'multiply_the_second_image_by')}\n",
         f"Measurement:{get_field(input_params,'image_or_measurement/measurement')}"
         ])
    return(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--inputs',
        help='JSON inputs from Galaxy'
    )

    parser.add_argument(
        '-p', '--pipeline',
        help='CellProfiler pipeline'
    )

    args = parser.parse_args()

    input_json_path = args.inputs
    input_pipeline_file = args.pipeline
    
    input_params_galaxy = json.load(open(input_json_path, "r"))

    pipeline_lines = get_pipeline_lines(input_pipeline_file)
    current_module_num = get_total_number_of_modules(pipeline_lines) + 1
    pipeline_lines = update_module_count(pipeline_lines, current_module_num)
    module_header = build_header(current_module_num, MODULE_NAME)
    main_block = build_main_block(input_params_galaxy)
    module_chunk = f"\n{module_header}{main_block}"
    if ("first_image" in input_params_galaxy['operation']):
        first_image_block = build_first_image_block(
            input_params_galaxy['operation']['first_image'])
        module_chunk = "\n".join([module_chunk, first_image_block])
    if ("second_image" in input_params_galaxy['operation']):
        second_image_block = build_second_image_block(
            input_params_galaxy['operation']['second_image'])
        module_chunk = "\n".join([module_chunk, second_image_block])

    pipeline_lines.append(module_chunk)

    write_pipeline(OUTPUT_FILENAME, pipeline_lines)
