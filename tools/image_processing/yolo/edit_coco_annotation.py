import argparse
import json
import sys
from copy import deepcopy
from pathlib import Path


def parse_ids(raw):
    ids = set()
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            a, b = token.split("-", 1)
            ids.update(range(int(a), int(b) + 1))
        else:
            ids.add(int(token))
    return ids


def parse_int_or_none(s):
    s = s.strip()
    return int(s) if s else None


def build_frame_map(images):
    def frame_key(img):
        return img.get("frame_index", img.get("frame_id", img["id"]))

    sorted_imgs = sorted(images, key=frame_key)
    return {
        img["id"]: img.get("frame_index", img.get("frame_id", i))
        for i, img in enumerate(sorted_imgs)
    }


def in_frame_range(ann, frame_map, frame_min, frame_max):
    frame = frame_map.get(ann["image_id"])
    if frame is None:
        return True
    if frame_min is not None and frame < frame_min:
        return False
    if frame_max is not None and frame > frame_max:
        return False
    return True


def parse_entries_keep(raw_entries):
    entries = []
    for raw in raw_entries:
        parts = raw.split("|")
        if len(parts) != 4:
            print(f"Invalid entry for keep mode: '{raw}'", file=sys.stderr)
            sys.exit(1)
        ids_str, rename, fmin_str, fmax_str = parts
        entries.append(
            {
                "ids": parse_ids(ids_str),
                "rename": rename.strip() or None,
                "frame_min": parse_int_or_none(fmin_str),
                "frame_max": parse_int_or_none(fmax_str),
            }
        )
    return entries


def parse_entries_remove(raw_entries):
    entries = []
    for raw in raw_entries:
        parts = raw.split("|")
        if len(parts) != 3:
            print(f"Invalid entry for remove mode: '{raw}'", file=sys.stderr)
            sys.exit(1)
        ids_str, fmin_str, fmax_str = parts
        entries.append(
            {
                "ids": parse_ids(ids_str),
                "frame_min": parse_int_or_none(fmin_str),
                "frame_max": parse_int_or_none(fmax_str),
            }
        )
    return entries


def parse_entries_rename(raw_entries):
    entries = []
    for raw in raw_entries:
        parts = raw.split("|")
        if len(parts) != 2:
            print(f"Invalid entry for rename mode: '{raw}'", file=sys.stderr)
            sys.exit(1)
        ids_str, name = parts
        entries.append(
            {
                "ids": parse_ids(ids_str),
                "name": name.strip(),
            }
        )
    return entries


def _apply_renames(
    annotations, categories, track_to_entry, all_track_ids, name_key="rename"
):
    name_to_cat_id = {cat["name"]: cat["id"] for cat in categories}
    max_cat_id = max((cat["id"] for cat in categories), default=0)
    track_to_new_cat = {}

    for entry in track_to_entry.values():
        new_name = entry.get(name_key)
        if not new_name:
            continue
        for tid in entry["ids"]:
            if tid not in all_track_ids:
                continue
            if new_name not in name_to_cat_id:
                max_cat_id += 1
                categories.append(
                    {"id": max_cat_id, "name": new_name, "supercategory": ""}
                )
                name_to_cat_id[new_name] = max_cat_id
                print(f"  New category: id={max_cat_id}, name='{new_name}'")
            else:
                print(
                    f"  Reusing category: id={name_to_cat_id[new_name]}, name='{new_name}'"
                )
            track_to_new_cat[tid] = name_to_cat_id[new_name]

    for ann in annotations:
        t = ann.get("track_id")
        if t is not None and t in track_to_new_cat:
            ann["category_id"] = track_to_new_cat[t]


def _cleanup(output, annotations, categories, images):
    used_cat_ids = {a["category_id"] for a in annotations}
    categories = [cat for cat in categories if cat["id"] in used_cat_ids]

    used_image_ids = {a["image_id"] for a in annotations}
    images = [img for img in images if img["id"] in used_image_ids]

    for new_id, ann in enumerate(annotations, start=1):
        ann["id"] = new_id

    output["annotations"] = annotations
    output["categories"] = categories
    output["images"] = images

    print(
        f"\nResult: {len(annotations)} annotation(s), {len(images)} image(s) kept."
    )
    return output


def _build_track_map(entries):
    """Last entry wins when an ID appears in multiple entries."""
    track_to_entry = {}
    for entry in entries:
        for tid in entry["ids"]:
            track_to_entry[tid] = entry
    return track_to_entry


def _print_frame_range(frame_min, frame_max):
    if frame_min is None and frame_max is None:
        return
    if frame_min is not None and frame_max is not None:
        print(
            f"  You chose to keep frames {frame_min} to {frame_max} (inclusive)"
        )
    elif frame_min is not None:
        print(f"  You chose to keep frames from {frame_min} onwards")
    else:
        print(f"  You chose to keep frames up to {frame_max} (inclusive)")


def edit_coco_keep(data, entries):
    output = deepcopy(data)
    images = output["images"]
    annotations = output["annotations"]
    categories = output["categories"]

    frame_map = build_frame_map(images)
    with_track = [a for a in annotations if a.get("track_id") is not None]
    without_track = [a for a in annotations if a.get("track_id") is None]

    all_track_ids = {a["track_id"] for a in with_track}
    track_to_entry = _build_track_map(entries)

    unlisted = all_track_ids - set(track_to_entry)
    if unlisted:
        print(f"  Removing unlisted track IDs: {sorted(unlisted)}")

    for entry in entries:
        _print_frame_range(entry["frame_min"], entry["frame_max"])

    kept = []
    for ann in with_track:
        tid = ann["track_id"]
        entry = track_to_entry.get(tid)
        print(f"tid : {tid}")
        print(f"entry : {entry}")
        if entry is None:
            continue
        frame_val = frame_map.get(ann["image_id"])
        print(
            f"  ann image_id={ann['image_id']} -> frame_index={frame_val} | frame_min={entry['frame_min']} frame_max={entry['frame_max']} -> kept={in_frame_range(ann, frame_map, entry['frame_min'], entry['frame_max'])}"
        )
        if not in_frame_range(
            ann, frame_map, entry["frame_min"], entry["frame_max"]
        ):
            continue
        kept.append(ann)

    if any(
        e["frame_min"] is not None or e["frame_max"] is not None
        for e in entries
    ):
        kept_frames = sorted(
            {
                frame_map[ann["image_id"]]
                for ann in kept
                if ann["image_id"] in frame_map
            }
        )
        print(f"  Frames in output: {kept_frames}")

    _apply_renames(
        kept, categories, track_to_entry, all_track_ids, name_key="rename"
    )
    return _cleanup(output, kept + without_track, categories, images)


def edit_coco_remove(data, entries):
    output = deepcopy(data)
    images = output["images"]
    annotations = output["annotations"]
    categories = output["categories"]

    frame_map = build_frame_map(images)
    with_track = [a for a in annotations if a.get("track_id") is not None]
    without_track = [a for a in annotations if a.get("track_id") is None]

    track_to_entry = _build_track_map(entries)

    for entry in entries:
        _print_frame_range(entry["frame_min"], entry["frame_max"])

    kept = []
    for ann in with_track:
        tid = ann["track_id"]
        entry = track_to_entry.get(tid)
        if entry is None:
            kept.append(ann)
            continue
        if not in_frame_range(
            ann, frame_map, entry["frame_min"], entry["frame_max"]
        ):
            kept.append(ann)

    removed = len(with_track) - len(kept)
    if removed:
        print(f"  {removed} annotation(s) removed.")

    return _cleanup(output, kept + without_track, categories, images)


def edit_coco_rename(data, entries):
    output = deepcopy(data)
    images = output["images"]
    annotations = output["annotations"]
    categories = output["categories"]

    all_track_ids = {
        a["track_id"] for a in annotations if a.get("track_id") is not None
    }
    track_to_entry = _build_track_map(entries)

    unknown = set(track_to_entry) - all_track_ids
    if unknown:
        print(f"  Warning: track ID(s) not found in file: {sorted(unknown)}")

    _apply_renames(
        annotations, categories, track_to_entry, all_track_ids, name_key="name"
    )
    return _cleanup(output, annotations, categories, images)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations", required=True)
    parser.add_argument(
        "--mode", required=True, choices=["keep", "remove", "rename"]
    )
    parser.add_argument("--entry", action="append", default=[])
    args = parser.parse_args()

    input_path = Path(args.annotations)
    if not input_path.exists():
        print(f"COCO file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    with open(input_path, encoding="utf-8") as f:
        data = json.load(f)

    print(f"Mode: {args.mode}")
    all_ids = {
        a["track_id"]
        for a in data.get("annotations", [])
        if a.get("track_id") is not None
    }
    print(f"Track IDs in file: {sorted(all_ids)}")

    if args.mode == "keep":
        entries = parse_entries_keep(args.entry)
        result = edit_coco_keep(data, entries)
    elif args.mode == "remove":
        entries = parse_entries_remove(args.entry)
        result = edit_coco_remove(data, entries)
    else:
        entries = parse_entries_rename(args.entry)
        result = edit_coco_rename(data, entries)

    output_path = Path("outputs") / "annotations_edited.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
