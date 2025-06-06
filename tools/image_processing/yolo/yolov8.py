import os
from ultralytics import YOLO
import shutil
import argparse
import numpy as np
import pathlib
from argparse import RawTextHelpFormatter
from termcolor import colored
import time
import cv2
from tifffile import imwrite
from collections import defaultdict


#
# Input arguments
#
parser = argparse.ArgumentParser(
    description='train/predict dataset with YOLOv8',
    epilog="""USAGE EXAMPLE:\n\n~~~~Prediction~~~~\n\
        python yolov8.py --test_path=/g/group/user/data --model_path=/g/cba/models --model_name=yolov8n --save_dir=/g/group/user/results --iou=0.7 --confidence=0.5 --image_size=320 --run_dir=/g/group/user/runs --foldername=batch --headless --num_classes=1 max_det=1 --class_names_file=/g/group/user/class_names.txt\n\
        \n~~~~Training~~~~ \n\
        python yolov8.py --train --yaml_path=/g/group/user/example.yaml  --model_path=/g/cba/models --model_name=yolov8n --run_dir=/g/group/user/runs/ --image_size=320 --epochs=150 --scale=0.3 --hsv_v=0.5 --model_format=pt --degrees=180 --class_names_file=/g/group/user/class_names.txt""", formatter_class=RawTextHelpFormatter)
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
parser.add_argument("--class_names_file",
                    help="Path to the text file containing class names.",
                    type=str)

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


#
# Functions
#
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

    train_save_path = os.path.expanduser('~/runs/' + args.mode + '/train/')
    if os.path.isdir(train_save_path):
        shutil.rmtree(train_save_path)
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
    # Remove prediction save path if already exists
    val_save_path = os.path.expanduser('~/runs/' + args.mode + '/val/')
    if os.path.isdir(val_save_path):
        shutil.rmtree(val_save_path)
    # Validate the model
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
        class_array = list(range(kwargs['num_classes']))
    else:
        class_array = [0, 1]

    if "max_det" in kwargs:
        maximum_detections = args.max_det
    else:
        maximum_detections = 300

    if "run_dir" in kwargs:
        run_save_dir = kwargs['run_dir']
    else:
        # Remove prediction save path if already exists
        pred_save_path = os.path.expanduser('~/runs/' + args.mode + '/predict/')
        if os.path.isdir(pred_save_path):
            shutil.rmtree(pred_save_path)
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
        train_save_path = os.path.expanduser('~/runs/' + args.mode + '/')
        if os.path.isfile(os.path.join(train_save_path,
                                       "train", "weights", "best.pt")) and (args.model_name == 'sam'):
            model = YOLO(os.path.join(train_save_path,
                                      "train", "weights", "best.pt"))
        else:
            model = YOLO(os.path.join(args.model_path,
                                      args.model_name + ".pt"))
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
        elif (args.mode == "track"):
            results = model.track(source=datapath_for_prediction,
                                  tracker=args.tracker_file,
                                  conf=args.confidence,
                                  iou=args.iou,
                                  persist=False,
                                  show=True,
                                  save=True,
                                  project=args.run_dir,
                                  name=args.foldername)
            # Store the track history
            track_history = defaultdict(lambda: [])

            for result in results:
                # Get the boxes and track IDs
                if result.boxes and result.boxes.is_track:
                    boxes = result.boxes.xywh.cpu()
                    track_ids = result.boxes.id.int().cpu().tolist()
                    # Visualize the result on the frame
                    frame = result.plot()
                    # Plot the tracks
                    for box, track_id in zip(boxes, track_ids):
                        x, y, w, h = box
                        track = track_history[track_id]
                        track.append((float(x), float(y)))  # x, y center point
                        if len(track) > 30:  # retain 30 tracks for 30 frames
                            track.pop(0)

                        # Draw the tracking lines
                        points = np.hstack(track).astype(np.int32).reshape((-1, 1, 2))
                        cv2.polylines(frame, [points], isClosed=False, color=(230, 230, 230), thickness=2)

                    # Display the annotated frame
                    cv2.imshow("YOLO11 Tracking", frame)
                    print(colored(f"Tracking results saved in : '{args.save_dir}' \n", 'green'))
        elif (args.mode == "segment"):
            # Read class names from the file
            with open(args.class_names_file, 'r') as f:
                class_names = [line.strip() for line in f.readlines()]
            # Create a mapping from class names to indices
            class_to_index = {class_name: i for i, class_name in enumerate(class_names)}

            # Save polygon coordinates
            for result in predictions:
                # Create binary mask
                img = np.copy(result.orig_img)
                filename = pathlib.Path(result.path).stem
                b_mask = np.zeros(img.shape[:2], np.uint8)
                mask_save_as = str(pathlib.Path(os.path.join(args.save_dir, filename + "_mask.tif")).absolute())
                # Define output file path for text file
                output_filename = os.path.splitext(filename)[0] + ".txt"
                txt_save_as = str(pathlib.Path(os.path.join(args.save_dir, filename + ".txt")).absolute())

                for c, ci in enumerate(result):
                    #  Extract contour result
                    contour = ci.masks.xy.pop()
                    #  Changing the type
                    contour = contour.astype(np.int32)
                    #  Reshaping
                    contour = contour.reshape(-1, 1, 2)
                    # Draw contour onto mask
                    _ = cv2.drawContours(b_mask, [contour], -1, (255, 255, 255), cv2.FILLED)

                    # Normalized polygon points
                    points = ci.masks.xyn.pop()
                    obj_class = int(ci.boxes.cls.to("cpu").numpy().item())
                    confidence = result.boxes.conf.to("cpu").numpy()[c]

                    with open(txt_save_as, 'a') as f:
                        segmentation_points = ['{} {}'.format(points[i][0], points[i][1]) for i in range(len(points))]
                        segmentation_points_string = ' '.join(segmentation_points)
                        line = '{} {} {}\n'.format(obj_class, segmentation_points_string, confidence)
                        f.write(line)

                imwrite(mask_save_as, b_mask, imagej=True)  # save image
                print(colored(f"Saved cropped image as : \n '{mask_save_as}' \n", 'magenta'))
                print(colored(f"Polygon coordinates saved as : \n '{txt_save_as}' \n", 'cyan'))

        else:
            raise Exception(("Currently only 'detect' and 'segment' modes are available"))
