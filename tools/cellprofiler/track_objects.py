#!/usr/bin/env python

import argparse
import json

from cp_common_functions import get_json_value
from cp_common_functions import get_pipeline_lines
from cp_common_functions import get_total_number_of_modules
from cp_common_functions import write_pipeline
from cp_common_functions import update_module_count
from cp_common_functions import INDENTATION

MODULE_NAME = "TrackObjects"
OUTPUT_FILENAME = "output.cppipe"


def build_header(module_name, module_number):
    result = "|".join([f"{module_name}:[module_num:{module_number}",
                       "svn_version:\\'Unknown\\'",
                       "variable_revision_number:7",
                       "show_window:True",
                       "notes:\\x5B\\'Track the embryos across images using the Overlap method\\x3A tracked objects are identified by the amount of frame-to-frame overlap. Save an image of embryos labeled with a unique number across time.\\'\\x5D",
                       "batch_state:array(\\x5B\\x5D, dtype=uint8)",
                       "enabled:True",
                       "wants_pause:False]\n"])
    return result


def build_main_block(input_params):
    result = INDENTATION.join([f"{INDENTATION}Choose a tracking method:{get_json_value(input_params,'con_tracking_method.tracking_method')}\n",
                               f"Select the objects to track:{get_json_value(input_params,'object_to_track')}\n"
                               ])

    tracking_method = get_json_value(input_params, 'con_tracking_method.tracking_method')

    obj_measurement = "None"  # default value
    if tracking_method == "Measurements":
        measurement_category = get_json_value(input_params, 'con_tracking_method.con_tracking_category.measurement_category')
        measurement = get_json_value(input_params, 'con_tracking_method.con_tracking_category.measurement')

        if measurement_category == "Intensity" or measurement_category == "Location":
            img_measure = get_json_value(input_params, 'con_tracking_method.con_tracking_category.img_measure')
            obj_measurement = f"{measurement_category}_{measurement}_{img_measure}"
        else:
            obj_measurement = f"{measurement_category}_{measurement}"

    result += INDENTATION.join([f"{INDENTATION}Select object measurement to use for tracking:{obj_measurement}\n"])

    if tracking_method == "LAP":  # no max distance required, set default for pipeline
        max_distance = 50
    else:
        max_distance = get_json_value(input_params, 'con_tracking_method.max_distance')

    result += INDENTATION.join([f"{INDENTATION}Maximum pixel distance to consider matches:{max_distance}\n"])

    display_option = get_json_value(input_params, 'con_tracking_method.display_option')

    output_img_name = "TrackedCells"  # default value, required by cppipe regardless of its presence in UI
    save = get_json_value(input_params, 'con_tracking_method.con_save_coded_img.save_coded_img')
    if save == "Yes":
        output_img_name = get_json_value(input_params, 'con_tracking_method.con_save_coded_img.name_output_img')

    result += INDENTATION.join(
        [f"{INDENTATION}Select display option:{display_option}\n",
         f"Save color-coded image?:{save}\n",
         f"Name the output image:{output_img_name}\n"
         ])

    # LAP method default values
    movement_model = "Both"
    no_std = 3.0
    radius_limit_max = 10.0
    radius_limit_min = 2.0
    radius = "2.0,10.0"
    run_second = "Yes"
    gap_closing = 40
    split_alt = 40
    merge_alt = 40
    max_gap_displacement = 5
    max_split = 50
    max_merge = 50
    max_temporal = 5
    max_mitosis_dist = 40
    mitosis_alt = 80

    # LAP method
    if tracking_method == "LAP":
        movement_model = get_json_value(input_params, 'con_tracking_method.movement_model')
        no_std = get_json_value(input_params, 'con_tracking_method.movement_model')
        radius_limit_max = get_json_value(input_params, 'con_tracking_method.max_radius')
        radius_limit_min = get_json_value(input_params, 'con_tracking_method.min_radius')
        radius = f"{radius_limit_min},{radius_limit_max}"

        run_second = get_json_value(input_params, 'con_tracking_method.con_second_lap.second_lap')
        if run_second == "Yes":
            gap_closing = get_json_value(input_params, 'con_tracking_method.con_second_lap.gap_closing')
            split_alt = get_json_value(input_params, 'con_tracking_method.con_second_lap.split_alt')
            merge_alt = get_json_value(input_params, 'con_tracking_method.con_second_lap.merge_alt')
            max_gap_displacement = get_json_value(input_params, 'con_tracking_method.con_second_lap.max_gap_displacement')
            max_split = get_json_value(input_params, 'con_tracking_method.con_second_lap.max_split')
            max_merge = get_json_value(input_params, 'con_tracking_method.con_second_lap.max_merge')
            max_temporal = get_json_value(input_params, 'con_tracking_method.con_second_lap.max_temporal')
            max_mitosis_dist = get_json_value(input_params, 'con_tracking_method.con_second_lap.max_mitosis_distance')
            mitosis_alt = get_json_value(input_params, 'con_tracking_method.con_second_lap.mitosis_alt')

    result += INDENTATION.join(
        [f"{INDENTATION}Select the movement model:{movement_model}\n",
         f"Number of standard deviations for search radius:{no_std}\n",
         f"Search radius limit, in pixel units (Min,Max):{radius}\n",
         f"Run the second phase of the LAP algorithm?:{run_second}\n",
         f"Gap closing cost:{gap_closing}\n",
         f"Split alternative cost:{split_alt}\n",
         f"Merge alternative cost:{merge_alt}\n",
         f"Maximum gap displacement, in pixel units:{max_gap_displacement}\n",
         f"Maximum split score:{max_split}\n",
         f"Maximum merge score:{max_merge}\n",
         f"Maximum temporal gap, in frames:{max_temporal}\n"
         ])

    # common section
    filter_by_lifetime = get_json_value(input_params, 'con_tracking_method.con_filter_by_lifetime.filter_by_lifetime')
    use_min = "Yes"  # default
    min_life = 1  # default
    use_max = "No"  # default
    max_life = 100  # default

    if filter_by_lifetime == "Yes":
        use_min = get_json_value(input_params, 'con_tracking_method.con_fileter_by_lifetime.con_use_min.use_min')
        if use_min == "Yes":
            min_life = get_json_value(input_params, 'con_tracking_method.con_fileter_by_lifetime.con_use_min.min_lifetime')

        use_max = get_json_value(input_params, 'con_tracking_method.con_fileter_by_lifetime.con_use_max.use_max')
        if use_max == "Yes":
            max_life = get_json_value(input_params, 'con_tracking_method.con_fileter_by_lifetime.con_use_max_lifetime')

    result += INDENTATION.join(
        [f"{INDENTATION}Filter objects by lifetime?:{filter_by_lifetime}\n",
         f"Filter using a minimum lifetime?:{use_min}\n",
         f"Minimum lifetime:{min_life}\n",
         f"Filter using a maximum lifetime?:{use_max}\n",
         f"Maximum lifetime:{max_life}\n"
         ])

    # print 2 leftover from LAP
    result += INDENTATION.join(
        [f"{INDENTATION}Mitosis alternative cost:{mitosis_alt}\n",
         f"Maximum mitosis distance, in pixel units:{max_mitosis_dist}\n"
         ])

    # Follow Neighbors
    # defaults
    avg_cell_diameter = 35.0
    use_adv = "No"
    cost_of_cell = 15.0
    weight_of_area_diff = 25.0

    if tracking_method == "Follow Neighbors":
        avg_cell_diameter = get_json_value(input_params, 'con_tracking_method.avg_diameter')
        use_adv = get_json_value(input_params, 'con_tracking_method.con_use_adv_parameter.use_adv_parameter')
        if use_adv == "Yes":
            cost_of_cell = get_json_value(input_params, 'con_tracking_method.con_use_adv_parameter.cost')
            weight_of_area_diff = get_json_value(input_params, 'con_tracking_method.con_use_adv_parameter.weight')

    result += INDENTATION.join(
        [f"{INDENTATION}Average cell diameter in pixels:{avg_cell_diameter}\n",
         f"Use advanced configuration parameters:{use_adv}\n",
         f"Cost of cell to empty matching:{cost_of_cell}\n",
         f"Weight of area difference in function matching cost:{weight_of_area_diff}"
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
