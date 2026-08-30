import argparse
import csv
import os
import pathlib
import time
from argparse import RawTextHelpFormatter
from collections import defaultdict

import cv2
import numpy as np
from termcolor import colored
from tifffile import imwrite
from ultralytics import YOLO

#
# Input arguments
#
parser = argparse.ArgumentParser(
    description='train/predict dataset with YOLOv8',
    epilog="""USAGE EXAMPLE:\n\n~~~~Prediction~~~~\n\
        python yolov8.py --test_path=/g/group/user/data --model_path=/g/cba/models --model_name=yolov8n --save_dir=/g/group/user/results --iou=0.7 --confidence=0.5 --image_size=320 --run_dir=/g/group/user/runs --foldername=batch --headless --num_classes=1 max_det=1 \n\
        \n~~~~Training~~~~ \n\
        python yolov8.py --train --yaml_path=/g/group/user/example.yaml  --model_path=/g/cba/models --model_name=yolov8n --run_dir=/g/group/user/runs/ --image_size=320 --epochs=150 --scale=0.3 --hsv_v=0.5 --model_format=pt --degrees=180 """, formatter_class=RawTextHelpFormatter)
parser.add_argument("--dir_path",
                    help=(
                        "Path to the training data directory."
                    ),
                    type=str)
parser.add_argument("--yaml_path",
                    help=(
                        "YAML file with all the data paths"
                        " i.e. for train, test, valid data."
                    ),
                    type=str)
parser.add_argument("--test_path",
                    help=(
                        "Path to the prediction folder."
                    ),
                    type=str)
parser.add_argument("--save_dir",
                    help=(
                        "Path to the directory where bounding boxes text files"
                        " would be saved."
                    ),
                    type=str)
parser.add_argument("--run_dir",
                    help=(
                        "Path where overlaid images would be saved."
                        "For example: `RUN_DIR=projectName/results`."
                        "This should exist."
                    ),
                    type=str)
parser.add_argument("--foldername",
                    help=("Folder to save overlaid images.\n"
                          "For example: FOLDERNAME=batch.\n"
                          "This should not exist as a new folder named `batch`\n"
                          " will be created in RUN_DIR.\n"
                          " If it exists already then, a new folder named `batch1`\n"
                          " will be created automatically as it does not overwrite\n"
                          ),
                    type=str)

# For selecting and loading model
parser.add_argument("--model_name",
                    help=("Models for task `detect` can be seen here:\n"
                          "https://docs.ultralytics.com/tasks/detect/#models \n\n"
                          "Models for task `segment` can be seen here:\n"
                          "https://docs.ultralytics.com/tasks/segment/#models \n\n"
                          " . Use `yolov8n` for `detect` tasks. "
                          "For custom model, use `best`"
                          ),
                    default='yolov8n', type=str)
parser.add_argument("--model_path",
                    help="Full absolute path to the model directory",
                    type=str)
parser.add_argument("--model_format",
                    help="Format of the YOLO model i.e pt, yaml etc.",
                    default='pt', type=str)
# For training the model and prediction
parser.add_argument("--mode",
                    help=(
                        "detection, segmentation, classification, and pose \n. "
                        " Only detection mode available currently i.e. `detect`"
                    ), default='detect', type=str)
parser.add_argument('--train',
                    help="Do training",
                    action='store_true')
parser.add_argument("--confidence",
                    help="Confidence value (0-1) for each detected bounding box",
                    default=0.5, type=float)
parser.add_argument("--epochs",
                    help="Number of epochs for training. Default: 100",
                    default=100, type=int)
parser.add_argument("--init_lr",
                    help="Number of epochs for training. Default: 100",
                    default=0.01, type=float)
parser.add_argument("--weight_decay",
                    help="Number of epochs for training. Default: 100",
                    default=0.0005, type=float)

parser.add_argument("--num_classes",
                    help="Number of classes to be predicted. Default: 2",
                    default=2, type=int)
parser.add_argument("--iou",
                    help="Intersection over union (IoU) threshold for NMS",
                    default=0.7, type=float)
parser.add_argument("--image_size",
                    help=("Size of input image to be used only as integer of w,h. \n"
                          "For training choose <= 1000. \n\n"
                          "Prediction will be done on original image size"
                          ),
                    default=320, type=int)
parser.add_argument("--max_det",
                    help=("Maximum number of detections allowed per image. \n"
                          "Limits the total number of objects the model can detect in a single inference, \n"
                          "preventing excessive outputs in dense scenes.\n\n"
                          ),
                    default=300, type=int)
# For tracking
parser.add_argument("--tracker_file",
                    help=("Path to the configuration file of the tracker used. \n"),
                    default='bytetrack.yaml', type=str)

# For headless operation
parser.add_argument('--headless', action='store_true')
parser.add_argument('--nextflow', action='store_true')


# For data augmentation
parser.add_argument("--hsv_h",
                    help="(float) image HSV-Hue augmentation (fraction)",
                    default=0.015, type=float)
parser.add_argument("--hsv_s",
                    help="(float) image HSV-Saturation augmentation (fraction)",
                    default=0.7, type=float)
parser.add_argument("--hsv_v",
                    help="(float) image HSV-Value augmentation (fraction)",
                    default=0.4, type=float)
parser.add_argument("--degrees",
                    help="(float) image rotation (+/- deg)",
                    default=0.0, type=float)
parser.add_argument("--translate",
                    help="(float) image translation (+/- fraction)",
                    default=0.1, type=float)
parser.add_argument("--scale",
                    help="(float) image scale (+/- gain)",
                    default=0.5, type=float)
parser.add_argument("--shear",
                    help="(float) image shear (+/- deg)",
                    default=0.0, type=float)
parser.add_argument("--perspective",
                    help="(float) image perspective (+/- fraction), range 0-0.001",
                    default=0.0, type=float)
parser.add_argument("--flipud",
                    help="(float) image flip up-down (probability)",
                    default=0.0, type=float)
parser.add_argument("--fliplr",
                    help="(float) image flip left-right (probability)",
                    default=0.5, type=float)
parser.add_argument("--mosaic",
                    help="(float) image mosaic (probability)",
                    default=1.0, type=float)
parser.add_argument("--crop_fraction",
                    help="(float) crops image to a fraction of its size to "
                    "emphasize central features and adapt to object scales, "
                    "reducing background distractions",
                    default=1.0, type=float)


# Train a new model on the dataset mentioned in yaml file
def trainModel(model_path, model_name, yaml_filepath, **kwargs):
    if "imgsz" in kwargs:
        image_size = kwargs['imgsz']
    else:
        image_size = 320

    if "epochs" in kwargs:
        n_epochs = kwargs['epochs']
    else:
        n_epochs = 100

    if "hsv_h" in kwargs:
        aug_hsv_h = kwargs['hsv_h']
    else:
        aug_hsv_h = 0.015

    if "hsv_s" in kwargs:
        aug_hsv_s = kwargs['hsv_s']
    else:
        aug_hsv_s = 0.7

    if "hsv_v" in kwargs:
        aug_hsv_v = kwargs['hsv_v']
    else:
        aug_hsv_v = 0.4

    if "degrees" in kwargs:
        aug_degrees = kwargs['degrees']
    else:
        aug_degrees = 10.0

    if "translate" in kwargs:
        aug_translate = kwargs['translate']
    else:
        aug_translate = 0.1

    if "scale" in kwargs:
        aug_scale = kwargs['scale']
    else:
        aug_scale = 0.2

    if "shear" in kwargs:
        aug_shear = kwargs['shear']
    else:
        aug_shear = 0.0

    if "shear" in kwargs:
        aug_shear = kwargs['shear']
    else:
        aug_shear = 0.0

    if "perspective" in kwargs:
        aug_perspective = kwargs['perspective']
    else:
        aug_perspective = 0.0

    if "fliplr" in kwargs:
        aug_fliplr = kwargs['fliplr']
    else:
        aug_fliplr = 0.5

    if "flipud" in kwargs:
        aug_flipud = kwargs['flipud']
    else:
        aug_flipud = 0.0

    if "mosaic" in kwargs:
        aug_mosaic = kwargs['mosaic']
    else:
        aug_mosaic = 1.0

    if "crop_fraction" in kwargs:
        aug_crop_fraction = kwargs['crop_fraction']
    else:
        aug_crop_fraction = 1.0

    if "weight_decay" in kwargs:
        weight_decay = kwargs['weight_decay']
    else:
        weight_decay = 1.0

    if "init_lr" in kwargs:
        init_lr = kwargs['init_lr']
    else:
        init_lr = 1.0

    # Load a pretrained YOLO model (recommended for training)
    if args.model_format == 'pt':
        model = YOLO(os.path.join(model_path, model_name + "." + args.model_format))
    else:
        model = YOLO(model_name + "." + args.model_format)
    model.train(data=yaml_filepath, epochs=n_epochs, project=args.run_dir,
                imgsz=image_size, verbose=True, hsv_h=aug_hsv_h,
                hsv_s=aug_hsv_s, hsv_v=aug_hsv_v, degrees=aug_degrees,
                translate=aug_translate, shear=aug_shear, scale=aug_scale,
                perspective=aug_perspective, fliplr=aug_fliplr,
                flipud=aug_flipud, mosaic=aug_mosaic, crop_fraction=aug_crop_fraction,
                weight_decay=weight_decay, lr0=init_lr)
    return model


# Validate the trained model
def validateModel(model):
    metrics = model.val()  # no args needed, dataset & settings remembered
    metrics.box.map    # map50-95
    metrics.box.map50  # map50
    metrics.box.map75  # map75
    metrics.box.maps   # a list contains map50-95 of each category


# Do predictions on images/videos using trained/loaded model
def predict(model, source_datapath, **kwargs):
    if "imgsz" in kwargs:
        image_size = kwargs['imgsz']
    else:
        image_size = 320

    if "conf" in kwargs:
        confidence = kwargs['conf']
    else:
        confidence = 0.5

    if "iou" in kwargs:
        iou_value = kwargs['iou']
    else:
        iou_value = 0.5

    if "num_classes" in kwargs:
        class_array = list(model.names.keys())
    else:
        class_array = [0, 1]

    if "max_det" in kwargs:
        maximum_detections = args.max_det
    else:
        maximum_detections = 300

    run_save_dir = kwargs['run_dir']  # For Galaxy, run_save_dir is always provided via xml wrapper
    if "foldername" in kwargs:
        save_folder_name = kwargs['foldername']

    # infer on a local image or directory containing images/videos
    prediction = model.predict(source=source_datapath, save=True, stream=True,
                               conf=confidence, imgsz=image_size,
                               save_conf=True, iou=iou_value, max_det=maximum_detections,
                               classes=class_array, save_txt=False,
                               project=run_save_dir, name=save_folder_name, verbose=True)
    return prediction


# Save bounding boxes
def save_yolo_bounding_boxes_to_txt(predictions, save_dir):
    """
    Function to save YOLO bounding boxes to text files.

    Parameters:
    - predictions: List of results from YOLO model inference.
    - save_dir: Directory where the text files will be saved.
    """
    for result in predictions:
        result = result.to("cpu").numpy()
        # Using bounding_boxes, confidence_scores, and class_num which are defined in the list
        bounding_boxes = result.boxes.xyxy  # Bounding boxes in xyxy format
        confidence_scores = result.boxes.conf  # Confidence scores
        class_nums = result.boxes.cls  # Class numbers

        # Create save directory if it doesn't exist
        save_path = pathlib.Path(save_dir).absolute()
        save_path.mkdir(parents=True, exist_ok=True)

        # Construct filename for the text file
        image_filename = pathlib.Path(result.path).stem
        text_filename = save_path / f"{image_filename}.txt"

        # Write bounding boxes info into the text file
        with open(text_filename, 'w') as f:
            for i in range(bounding_boxes.shape[0]):
                x1, y1, x2, y2 = bounding_boxes[i]
                confidence = confidence_scores[i]
                class_num = int(class_nums[i])
                f.write(f'{class_num:01} {x1:06.2f} {y1:06.2f} {x2:06.2f} {y2:06.2f} {confidence:0.02} \n')
            print(colored(f"Bounding boxes saved in: {text_filename}", 'green'))


# Main code
if __name__ == '__main__':
    args = parser.parse_args()

    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    # Train/load model
    if (args.train):
        model = trainModel(args.model_path, args.model_name, args.yaml_path,
                           imgsz=args.image_size, epochs=args.epochs,
                           hsv_h=args.hsv_h, hsv_s=args.hsv_s, hsv_v=args.hsv_v,
                           degrees=args.degrees, translate=args.translate,
                           shear=args.shear, scale=args.scale,
                           perspective=args.perspective, fliplr=args.fliplr,
                           flipud=args.flipud, mosaic=args.mosaic)
        validateModel(model)
    else:
        t = time.time()
        model = YOLO(os.path.join(args.model_path, args.model_name + ".pt"))
        model.info(verbose=True)
        elapsed = time.time() - t
        print(colored(f"\nYOLO model loaded in : '{elapsed}' sec \n", 'white', 'on_yellow'))

    if (args.save_dir):
        # Do predictions (optionally show image results with bounding boxes)
        t = time.time()
        datapath_for_prediction = args.test_path
        # Extracting class names from the model
        class_names = model.names
        predictions = predict(model, datapath_for_prediction,
                              imgsz=args.image_size, conf=args.confidence,
                              iou=args.iou, run_dir=args.run_dir,
                              foldername=args.foldername, num_classes=args.num_classes, max_det=args.max_det)
        elapsed = time.time() - t
        print(colored(f"\nYOLO prediction done in : '{elapsed}' sec \n", 'white', 'on_cyan'))

        if (args.mode == "detect"):
            # Save bounding boxes
            save_yolo_bounding_boxes_to_txt(predictions, args.save_dir)

            # Loop over each result
            for result in predictions:
                img = np.copy(result.orig_img)
                image_filename = pathlib.Path(result.path).stem
                overlay_path = os.path.join(args.save_dir, f"{image_filename}_overlay.jpg")

                for box, cls, conf in zip(result.boxes.xyxy, result.boxes.cls, result.boxes.conf):
                    x1, y1, x2, y2 = map(int, box.tolist())
                    class_num = int(cls.item())
                    confidence = conf.item()
                    label = f"{class_names[class_num]} {confidence:.2f}" if class_names else f"{class_num} {confidence:.2f}"

                    cv2.rectangle(img, (x1, y1), (x2, y2), color=(0, 255, 0), thickness=2)
                    cv2.putText(img, label, (x1, y1 - 10), cv2.FONT_HERSHEY_SIMPLEX,
                                0.5, (0, 255, 0), thickness=1)

                cv2.imwrite(overlay_path, img)
                print(colored(f"Overlay image saved at: {overlay_path}", 'cyan'))
        elif (args.mode == "track"):
            results = model.track(source=datapath_for_prediction,
                                  tracker=args.tracker_file,
                                  conf=args.confidence,
                                  iou=args.iou,
                                  persist=True,
                                  show=False,
                                  save=True,
                                  project=args.run_dir,
                                  name=args.foldername)
            # Store the track history
            track_history = defaultdict(lambda: [])

            tsv_path = os.path.join(args.save_dir, "tracks.tsv")
            with open(tsv_path, "w", newline="") as tsvfile:
                writer = csv.writer(tsvfile, delimiter='\t')
                writer.writerow(['track_id', 'frame', 'class', 'centroid_x', 'centroid_y'])
                frame_idx = 0
                for result in results:
                    # Get the boxes and track IDs
                    if result.boxes and result.boxes.is_track:
                        track_ids = result.boxes.id.int().cpu().tolist()
                        labels = result.boxes.cls.int().cpu().tolist() if hasattr(result.boxes, "cls") else [0] * len(track_ids)
                        # Prepare mask image
                        img_shape = result.orig_shape if hasattr(result, "orig_shape") else result.orig_img.shape
                        mask = np.zeros(img_shape[:2], dtype=np.uint16)
                        # Check if polygons (masks) are available
                        if hasattr(result, "masks") and result.masks is not None and hasattr(result.masks, "xy"):
                            polygons = result.masks.xy
                            for i, (track_id, label) in enumerate(zip(track_ids, labels)):
                                if i < len(polygons):
                                    contour = polygons[i].astype(np.int32)
                                    contour = contour.reshape(-1, 1, 2)
                                    cv2.drawContours(mask, [contour], -1, int(track_id), cv2.FILLED)
                                    # Calculate centroid of the polygon
                                    M = cv2.moments(contour)
                                    if M["m00"] != 0:
                                        cx = float(M["m10"] / M["m00"])
                                        cy = float(M["m01"] / M["m00"])
                                    else:
                                        cx, cy = 0.0, 0.0
                                    writer.writerow([track_id, frame_idx, label, cx, cy])
                        else:
                            # Fallback to bounding boxes if polygons are not available
                            boxes = result.boxes.xywh.cpu()
                            xyxy_boxes = result.boxes.xyxy.cpu().numpy()
                            for i, (box, xyxy, track_id, label) in enumerate(zip(boxes, xyxy_boxes, track_ids, labels)):
                                x, y, w, h = box
                                writer.writerow([track_id, frame_idx, label, float(x), float(y)])
                                x1, y1, x2, y2 = map(int, xyxy)
                                cv2.rectangle(mask, (x1, y1), (x2, y2), int(track_id), thickness=-1)
                        # Collect masks for TYX stack
                        if frame_idx == 0:
                            mask_stack = []
                        mask_stack.append(mask)
                    frame_idx += 1
            # Save TYX stack (T=frames, Y, X)
            if 'mask_stack' in locals() and len(mask_stack) > 0:
                tyx_array = np.stack(mask_stack, axis=0)
                # Remove string from last underscore in filename
                stem = pathlib.Path(result.path).stem
                stem = stem.rsplit('_', 1)[0] if '_' in stem else stem
                mask_save_as = str(pathlib.Path(os.path.join(args.save_dir, stem + "_mask.tiff")).absolute())
                imwrite(mask_save_as, tyx_array)
                print(colored(f"TYX mask stack saved as : '{mask_save_as}'", 'magenta'))
            print(colored(f"Tracking results saved in : '{args.save_dir}' \n", 'green'))
        elif (args.mode == "segment"):
            # Save polygon coordinates
            for result in predictions:
                # Create binary mask
                img = np.copy(result.orig_img)
                filename = pathlib.Path(result.path).stem
                b_mask = np.zeros(img.shape[:2], np.uint8)
                mask_save_as = str(pathlib.Path(os.path.join(args.save_dir, filename + "_mask.tiff")).absolute())
                # Define output file path for text file
                output_filename = os.path.splitext(filename)[0] + ".txt"
                txt_save_as = str(pathlib.Path(os.path.join(args.save_dir, filename + ".txt")).absolute())
                instance_id = 1  # Start instance IDs from 1
                for c, ci in enumerate(result):
                    # Extract contour result
                    contour = ci.masks.xy.pop()
                    contour = contour.astype(np.int32)
                    contour = contour.reshape(-1, 1, 2)
                    # Draw contour onto mask with unique instance id
                    _ = cv2.drawContours(b_mask, [contour], -1, instance_id, cv2.FILLED)

                    # Normalized polygon points
                    points = ci.masks.xyn.pop()
                    confidence = result.boxes.conf.to("cpu").numpy()[c]

                    with open(txt_save_as, 'a') as f:
                        segmentation_points = ['{} {}'.format(points[i][0], points[i][1]) for i in range(len(points))]
                        segmentation_points_string = ' '.join(segmentation_points)
                        line = '{} {} {}\n'.format(instance_id, segmentation_points_string, confidence)
                        f.write(line)

                    instance_id += 1  # Increment for next object

                imwrite(mask_save_as, b_mask, imagej=True)  # save label mask image
                print(colored(f"Saved label mask as : \n '{mask_save_as}' \n", 'magenta'))
                print(colored(f"Polygon coordinates saved as : \n '{txt_save_as}' \n", 'cyan'))
        else:
            raise Exception(("Currently only 'detect', 'segment' and 'track' modes are available"))
