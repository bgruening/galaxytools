import os
import shutil
import argparse
import json
from sklearn.model_selection import train_test_split

def get_basename(f):
    return os.path.splitext(os.path.basename(f))[0]

def pair_files(images_dir, labels_dir):

    img_files = [f for f in os.listdir(images_dir) ]
    lbl_files = [f for f in os.listdir(labels_dir) ]

    image_dict = {get_basename(f): f for f in img_files}
    label_dict = {get_basename(f): f for f in lbl_files}

    keys = sorted(set(image_dict) & set(label_dict))

    return [(image_dict[k], label_dict[k]) for k in keys]

def copy_pairs(pairs, image_src, label_src, image_dst, label_dst):
    os.makedirs(image_dst, exist_ok=True)
    os.makedirs(label_dst, exist_ok=True)
    for img, lbl in pairs:
        shutil.copy(os.path.join(image_src, img), os.path.join(image_dst, img))
        shutil.copy(os.path.join(label_src, lbl), os.path.join(label_dst, lbl))

def write_yolo_yaml(output_dir,meta_json):
    with open(meta_json, 'r') as f:
        meta = json.load(f)

    yolo_yaml_path = os.path.join(output_dir, "yolo.yml")
    with open(yolo_yaml_path, 'w') as f:
        f.write(f"path: {output_dir}\n")
        f.write(f"train: train\n")
        f.write(f"val: valid\n")
        f.write(f"test: test\n")
        f.write(f"\n")
        f.write(f"nc: {meta['training_params']['num_class']}\n")
        f.write(f"names: ['{meta['ds_name']}']\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--images", required=True)
    parser.add_argument("-y","--labels", required=True)
    parser.add_argument("-o","--output", required=True)
    parser.add_argument("-p","--train_percent", type=int, default=70)
    parser.add_argument("-m","--meta", required=True)
    args = parser.parse_args()

    all_pairs = pair_files(args.images, args.labels)
    train_size = args.train_percent / 100.0
    val_test_size = 1.0 - train_size

    train_pairs, val_test_pairs = train_test_split(all_pairs, test_size=val_test_size, random_state=42)
    val_pairs, test_pairs = train_test_split(val_test_pairs, test_size=0.5, random_state=42)

    copy_pairs(train_pairs, args.images, args.labels, os.path.join(args.output, "train/images"), os.path.join(args.output, "train/labels"))
    copy_pairs(val_pairs, args.images, args.labels, os.path.join(args.output, "valid/images"), os.path.join(args.output, "valid/labels"))
    copy_pairs(test_pairs, args.images, args.labels, os.path.join(args.output, "test/images"), os.path.join(args.output, "test/labels"))

    write_yolo_yaml(args.output,args.meta)

if __name__ == "__main__":
    main()
