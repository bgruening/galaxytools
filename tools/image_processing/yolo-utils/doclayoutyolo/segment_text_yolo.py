"""
Segment text using DocLayout Yolo model
"""

import argparse
import json
import os

import cv2
from doclayout_yolo import YOLOv10


def _as_numpy(value):
    if value is None:
        return None
    if hasattr(value, "cpu"):
        value = value.cpu()
    if hasattr(value, "numpy"):
        value = value.numpy()
    return value


def _class_name(names, class_id):
    if isinstance(names, dict):
        return names.get(class_id, names.get(str(class_id), str(class_id)))
    try:
        return names[class_id]
    except Exception:
        return str(class_id)


def _build_feature(x1, y1, x2, y2, properties):
    return {
        "type": "Feature",
        "geometry": {
            "type": "Polygon",
            "coordinates": [
                [
                    [x1, y1],
                    [x2, y1],
                    [x2, y2],
                    [x1, y2],
                    [x1, y1],
                ]
            ],
        },
        "properties": properties,
    }


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


def build_geojson(results):
    boxes = getattr(results, "boxes", None)
    if boxes is None or not len(boxes):
        return {"type": "FeatureCollection", "features": []}

    xyxy = _as_numpy(getattr(boxes, "xyxy", None))
    classes = _as_numpy(getattr(boxes, "cls", None))
    confidences = _as_numpy(getattr(boxes, "conf", None))
    names = getattr(results, "names", {})

    features = []
    for idx, bx in enumerate(xyxy):
        x1, y1, x2, y2 = (float(coord) for coord in bx[:4])
        properties = {"id": idx}

        if classes is not None and idx < len(classes):
            class_id = int(classes[idx])
            properties["class_id"] = class_id
            properties["class_name"] = _class_name(names, class_id)

        if confidences is not None and idx < len(confidences):
            properties["confidence"] = float(confidences[idx])

        features.append(_build_feature(x1, y1, x2, y2, properties))

    return {"type": "FeatureCollection", "features": features}


def extract_bb_crop(results, output_segmentation_coordinates):
    geojson_obj = build_geojson(results)

    with open(output_segmentation_coordinates, "w") as f:
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
        "-ie", "--input_image_ext", required=True, help="Input image file extension"
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
    args = vars(arg_parser.parse_args())
    model_path = args["yolo_model"]
    input_image_path = args["input_image"]
    input_ext = args["input_image_ext"]
    confidence = args["input_confidence"]
    image_size = args["input_image_size"]
    output_image_path = args["output_image"]
    output_segmentation_coordinates = args["output_geojson"]

    model_link = "yolo_model.pt"
    input_image = f"input_image.{input_ext}"
    output_image = f"output_image.{input_ext}"

    os.symlink(model_path, model_link)
    os.symlink(input_image_path, input_image)
    os.symlink(output_image_path, output_image)

    segmented_image = load_model_and_predict(
        model_link, input_image, confidence, image_size, output_image
    )
    extract_bb_crop(segmented_image, output_segmentation_coordinates)
