import os
import argparse
from skimage import exposure, img_as_ubyte
from skimage.io import imread, imsave

def apply_clahe(image):
    return exposure.equalize_adapthist(image, clip_limit=0.01)

def count_total_tif_files(folder_path):
    total_files = sum(
        len([f for f in files if f.lower().endswith('.tif')])
        for _, _, files in os.walk(folder_path) if files
    )
    return total_files

def process_images(input_folder, results_folder,clahe,convert_8bit):
    # Create the main "results" directory
    results_folder = os.path.join(results_folder, 'jpegs')
    os.makedirs(results_folder, exist_ok=True)

    # Count the total number of .tif files for the progress bar display
    total_files = count_total_tif_files(input_folder)
    # total_files = sum(len(files) for _, _, files in os.walk(input_folder) if files and any(f.endswith('.jpg') for f in files))

    for root, _, files in os.walk(input_folder):
        for file in files:
            if file.endswith('.tif'):
                # Derive the relative path to maintain the directory structure
                relative_path = os.path.relpath(root, input_folder)
                target_folder = os.path.join(results_folder, relative_path)
                os.makedirs(target_folder, exist_ok=True)

                file_path = os.path.join(root, file)
                image = imread(file_path)

                # Apply CLAHE
                if clahe:
                    enhanced_image = apply_clahe(image)

                # Convert to 8-bit
                if convert_8bit:
                    enhanced_image_8bit = img_as_ubyte(enhanced_image)

                # # Save as JPG
                jpg_filename = os.path.splitext(file)[0] + '.jpg'
                output_path = os.path.join(target_folder, jpg_filename)
                imsave(output_path, enhanced_image_8bit)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process .tif images with CLAHE and save as .jpg.",
                epilog="Usage example: python adjust_save_images.py input_folder=/g/data/ output_folder=/g/results")
    parser.add_argument("--input_folder", 
                        type=str, 
                        help="Path to the input directory containing .tif files.")
    parser.add_argument("--results_folder", 
                        type=str, 
                        help="Path to the output directory where processed images will be saved.")

    parser.add_argument("--clahe", action="store_true", help="Apply CLAHE to the images.")
    parser.add_argument("--convert8bit", action="store_true", help="Convert images to 8-bit before saving.")

    args = parser.parse_args()

    process_images(args.input_folder, args.results_folder,args.clahe, args.convert8bit)

