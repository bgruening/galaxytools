import argparse
import json
import os


def convert_json_to_yolo(input_dir, save_dir, class_names_file):

    with open(class_names_file, 'r') as f:
        class_names = [line.strip() for line in f.readlines()]

    class_to_index = {class_name: i for i, class_name in enumerate(class_names)}

    for filename in os.listdir(input_dir):
        json_filepath = os.path.join(input_dir, filename)

        with open(json_filepath, 'r') as f:
            data = json.load(f)

        image_width = data.get('imageWidth')
        image_height = data.get('imageHeight')
        if image_width is None or image_height is None:
            print(f"Skipping {filename}: missing image dimensions.")
            return

        annotations = data.get('shapes', [])

        base, _ = os.path.splitext(filename)
        output_file = f"{base}.txt"
        output_filepath = os.path.join(save_dir, output_file)

        with open(output_filepath, 'w') as f:
            for annotation in annotations:
                label = annotation.get('label')
                class_index = class_to_index[label]

                points = annotation.get('points', [])

                if not points:
                    print(f"No points found for annotation '{label}', skipping.")
                    continue

                x = [point[0] / image_width for point in points]
                y = [point[1] / image_height for point in points]

                segmentation_points = ['{} {}'.format(x[i], y[i]) for i in range(len(x))]
                segmentation_points_string = ' '.join(segmentation_points)
                line = '{} {}\n'.format(class_index, segmentation_points_string)
                f.write(line)

        print(f"Converted annotations saved to: {output_filepath}")


def main():
    parser = argparse.ArgumentParser(description="Convert JSON annotations to YOLO segment format.")
    parser.add_argument('-i', '--input_dir', type=str, help='Full path of the folder containing AnyLabeling JSON files.')
    parser.add_argument('-o', '--save_dir', type=str, help='Path to the directory to save converted YOLO files.')
    parser.add_argument('-c', '--class_names_file', type=str, help='Path to the text file containing class names, one per line.')
    args = parser.parse_args()

    convert_json_to_yolo(args.input_dir, args.save_dir, args.class_names_file)


if __name__ == "__main__":
    main()
