"""
Segment text using DocLayout Yolo model
"""

import argparse
import os

import cv2
import json
from doclayout_yolo import YOLOv10
from geojson import Feature, FeatureCollection
from shapely.geometry import box, mapping


def load_model_and_predict(
    model_path, input_image_path, input_confidence, image_size, output_image_path
):

    model = YOLOv10(model=model_path)

    det_res = model.predict(
        input_image_path, imgsz=int(image_size), conf=float(input_confidence)
    )
    annotated_frame = det_res[0].plot(pil=True, line_width=5, font_size=20)
    cv2.imwrite(output_image_path, annotated_frame)
    return det_res[0]


def extract_bb_crop(results, output_segmentation_coordiates):
    bounding_boxes = []
    features = []
    for bx in results.boxes.xyxy.cpu().numpy():
        x1, y1, x2, y2 = bx
        bounding_boxes.append((x1, y1, x2, y2))

    for i, (x1, y1, x2, y2) in enumerate(bounding_boxes):
        poly = box(x1, y1, x2, y2)
        feature = Feature(geometry=mapping(poly), properties={"id": i})
        features.append(feature)

    geojson_obj = FeatureCollection(features)

    # Save to file
    with open(output_segmentation_coordiates, "w") as f:
        json.dump(geojson_obj, f)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "-im", "--yolo_model", required=True, help="Input Yolo model"
    )
    arg_parser.add_argument(
        "-ii", "--input_image", required=True, help="Input image file"
    )
    arg_parser.add_argument(
        "-ic", "--input_confidence", required=True, help="Input confidence"
    )
    arg_parser.add_argument(
        "-is", "--input_image_size", required=True, help="Input image size"
    )
    arg_parser.add_argument("-oi", "--output_image", required=True, help="Output image")
    arg_parser.add_argument(
        "-ogj", "--output_geojson", required=True, help="Output segmented coordinates"
    )
    # get argument values
    args = vars(arg_parser.parse_args())
    model_path = args["yolo_model"]
    input_image_path = args["input_image"]
    confidence = args["input_confidence"]
    image_size = args["input_image_size"]
    output_image_path = args["output_image"]
    output_segmentation_coordiates = args["output_geojson"]

    model_link = "yolo_model.pt"
    input_image = "input_image.png"
    output_image = "output_image.png"

    os.symlink(model_path, model_link)
    os.symlink(input_image_path, input_image)
    os.symlink(output_image_path, output_image)

    segmented_image = load_model_and_predict(
        model_link, input_image, confidence, image_size, output_image
    )
    extract_bb_crop(segmented_image, output_segmentation_coordiates)
