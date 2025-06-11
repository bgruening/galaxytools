import argparse
import os

from skimage import exposure, img_as_ubyte
from skimage.io import imread, imsave


def apply_clahe(image):
    return exposure.equalize_adapthist(image, clip_limit=0.01)


def process_images(input_folder, results_folder):
    # Create the main "results" directory
    results_folder = os.path.join(results_folder, 'jpegs')
    os.makedirs(results_folder, exist_ok=True)

    for root, _, files in os.walk(input_folder):
        for file in files:
            if file.endswith(('.tif', '.tiff')):
                # Derive the relative path to maintain the directory structure
                relative_path = os.path.relpath(root, input_folder)
                target_folder = os.path.join(results_folder, relative_path)
                os.makedirs(target_folder, exist_ok=True)

                file_path = os.path.join(root, file)
                image = imread(file_path)
                enhanced_image = img_as_ubyte(image)

                # # Save as JPG
                jpg_filename = os.path.splitext(file)[0] + '.jpg'
                output_path = os.path.join(target_folder, jpg_filename)
                imsave(output_path, enhanced_image)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process .tiff images with CLAHE and/or save as .jpg.")
    parser.add_argument("--input_folder",
                        type=str,
                        help="Path to the input directory containing .tiff files.")
    parser.add_argument("--results_folder",
                        type=str,
                        help="Path to the output directory where processed images will be saved.")

    args = parser.parse_args()

    process_images(args.input_folder, args.results_folder)
