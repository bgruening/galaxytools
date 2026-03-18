#!/usr/bin/python3
import argparse
import copy
import json
import os
import os.path as osp

from tqdm import tqdm

DEFAULT_LABELME = {
    "version": "0.4.30",
    "flags": {},
    "shapes": [],
    "imagePath": None,
    "imageData": None,
    "imageHeight": None,
    "imageWidth": None,
    "text": ""
}

DEFAULT_SHAPE = {
    "label": None,
    "text": "",
    "points": [],
    "group_id": None,
    "shape_type": None,
    "flags": {}
}


def find_all_img_anns(coco):
    img_id_list = []
    anns_list = []
    for img_info in coco['images']:
        img_id_list.append(img_info['id'])
        anns_list.append([])
    for ann in tqdm(coco['annotations']):
        index = img_id_list.index(ann['image_id'])
        anns_list[index].append(ann)
    return coco['images'], anns_list


def coco2labelme(coco_path, outputs, anylabeling, custom_path):
    os.makedirs(outputs, exist_ok=True)
    with open(coco_path, 'r') as f:
        coco = json.loads(f.read())

    categories = {cat["id"]: cat["name"] for cat in coco["categories"]}
    img_info_list, anns_list = find_all_img_anns(coco)

    for img_info, anns in zip(img_info_list, anns_list):
        labelme_json = copy.deepcopy(DEFAULT_LABELME)
        labelme_json['imageHeight'] = img_info['height']
        labelme_json['imageWidth'] = img_info['width']

        shapes = []
        for ann in anns:
            if not ann.get('segmentation') or len(ann['segmentation']) == 0:
                continue

            shape = copy.deepcopy(DEFAULT_SHAPE)
            points = ann['segmentation'][0]
            points = [[float(x), float(y)] for x, y in zip(points[::2], points[1::2])]

            if len(points) == 1:
                shape_type = 'point'
            elif len(points) == 2:
                shape_type = 'line'
            elif len(points) < 1:
                continue
            else:
                shape_type = 'polygon'

            shape['points'] = points
            shape['shape_type'] = shape_type
            shape['label'] = categories[ann['category_id']]
            shapes.append(shape)

        labelme_json['shapes'] = shapes

        filename = osp.basename(img_info['file_name'])
        if anylabeling == "anylabeling":
            labelme_json['imagePath'] = osp.join("..", "input_images", filename)
        elif anylabeling == "simple":
            labelme_json['imagePath'] = osp.join("..", filename)
        else:
            labelme_json['imagePath'] = osp.join(custom_path, filename)

        print(f"labelme_json['imagePath'] {labelme_json['imagePath']}")

        out_filename = osp.splitext(osp.basename(img_info['file_name']))[0] + '.json'
        with open(osp.join(outputs, out_filename), 'w') as f:
            f.write(json.dumps(labelme_json, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--coco', type=str)
    parser.add_argument('--outputs', type=str)
    parser.add_argument('--anylabeling', type=str)
    parser.add_argument('--custom_path', type=str)
    opt = parser.parse_args()

    coco_files = opt.coco.split(",")
    for coco in coco_files:
        coco2labelme(coco, opt.outputs, opt.anylabeling, opt.custom_path)
