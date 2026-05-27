import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import cv2

import numpy as np

import matplotlib.colors as mcolors

from pycocotools import mask as mask_utils


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize COCO annotations on images"
    )
    parser.add_argument(
        "--annotations",
        required=True,
        help="Path to COCO JSON annotation file",
    )
    parser.add_argument(
        "--outdir",
        default="outputs",
        help="Output directory"
    )
    parser.add_argument(
        "--filter_categories",
        default="",
        help="Comma-separated category names to display (empty = all)",
    )
    parser.add_argument("--show_bbox", default="true")
    parser.add_argument("--show_mask", default="true")
    parser.add_argument("--show_labels", default="true")
    parser.add_argument(
        "--show_count",
        default="false",
        help="Draw annotation count in image corner",
    )
    parser.add_argument("--mask_opacity", type=float, default=0.4)
    parser.add_argument("--bbox_thickness", type=int, default=2)
    parser.add_argument(
        "--color_mode",
        default="per_category",
        choices=["per_category", "per_instance"],
    )
    parser.add_argument(
        "--output_format", default="png", choices=["png", "jpg"]
    )
    parser.add_argument("--font_scale", type=float, default=0.6)
    parser.add_argument(
        "--output_mode",
        default="frames",
        choices=["frames", "video", "both"],
        help="Output mode: frames (individual images), video (compiled video)"
        ", or both",
    )
    parser.add_argument(
        "--video_fps",
        type=float,
        default=25.0,
        help="Frame rate for video output",
    )
    return parser.parse_args()


def _bool(value: str) -> bool:
    return value.lower() == "true"


_TABLEAU_PALETTE: List[Tuple[int, int, int]] = [
    (int(b * 255), int(g * 255), int(r * 255))
    for hex_color in mcolors.TABLEAU_COLORS.values()
    for r, g, b in [mcolors.to_rgb(hex_color)]
]


def color_for_id(id: int) -> Tuple[int, int, int]:
    return _TABLEAU_PALETTE[id % len(_TABLEAU_PALETTE)]


def decode_segmentation(seg, height: int, width: int) -> Optional[np.ndarray]:
    """Decode a COCO polygon or RLE segmentation into a binary uint8 mask."""
    if not seg:
        return None
    if isinstance(seg, list):
        mask = np.zeros((height, width), dtype=np.uint8)
        for poly in seg:
            if len(poly) < 6:
                continue
            pts = (
                np.array(poly, dtype=np.float32)
                .reshape(-1, 2)
                .astype(np.int32)
            )
            cv2.fillPoly(mask, [pts], 1)
        return mask
    if isinstance(seg, dict):
        rle = mask_utils.frPyObjects(seg, height, width)
        return mask_utils.decode(rle).astype(np.uint8)
    return None


def draw_label(
    img: np.ndarray,
    text: str,
    x: int,
    y: int,
    color: Tuple[int, int, int],
    font_scale: float,
) -> None:
    thickness = 1
    (tw, th), baseline = cv2.getTextSize(
        text, cv2.FONT_HERSHEY_SIMPLEX, font_scale, thickness
    )
    pad = 3
    top = max(y - th - baseline - pad, 0)
    cv2.rectangle(img, (x, top), (x + tw + pad * 2, y + baseline), color, -1)
    cv2.putText(
        img,
        text,
        (x + pad, y),
        cv2.FONT_HERSHEY_SIMPLEX,
        font_scale,
        (255, 255, 255),
        thickness,
        cv2.LINE_AA,
    )


def draw_annotations(
    image: np.ndarray,
    annotations: List[dict],
    category_map: Dict[int, str],
    *,
    show_bbox: bool,
    show_mask: bool,
    show_labels: bool,
    mask_opacity: float,
    bbox_thickness: int,
    color_mode: str,
    font_scale: float,
    category_colors: Dict[int, Tuple[int, int, int]],
) -> np.ndarray:
    height, width = image.shape[:2]
    result = image.copy()

    for ann in annotations:
        cat_id = ann.get("category_id", 0)
        if color_mode == "per_category":
            color = category_colors.get(cat_id, (0, 200, 0))
        else:  # per_instance: prefer track_id, fallback to annotation id
            color = color_for_id(ann.get("track_id", ann.get("id", 0)))

        # Mask
        if show_mask and ann.get("segmentation"):
            mask = decode_segmentation(ann["segmentation"], height, width)
            if mask is not None:
                px = mask == 1
                result[px] = (
                    result[px] * (1 - mask_opacity)
                    + np.array(color) * mask_opacity
                ).astype(np.uint8)

        # Bounding box
        if show_bbox and "bbox" in ann:
            x, y, w, h = (int(v) for v in ann["bbox"])
            cv2.rectangle(
                result, (x, y), (x + w, y + h), color, bbox_thickness
            )

        # Label
        if show_labels and "bbox" in ann:
            label = category_map.get(cat_id, str(cat_id))
            if "track_id" in ann:
                label = f"{label} [{ann['track_id']}]"
            x, y = int(ann["bbox"][0]), int(ann["bbox"][1])
            draw_label(result, label, x, max(y - 2, 12), color, font_scale)

    return result


def draw_count_overlay(
    img: np.ndarray, count: int, font_scale: float
) -> np.ndarray:
    text = f"Annotations: {count}"
    (tw, th), baseline = cv2.getTextSize(
        text, cv2.FONT_HERSHEY_SIMPLEX, font_scale, 1
    )
    pad = 6
    cv2.rectangle(
        img, (0, 0), (tw + pad * 2, th + baseline + pad * 2), (0, 0, 0), -1
    )
    cv2.putText(
        img,
        text,
        (pad, th + pad),
        cv2.FONT_HERSHEY_SIMPLEX,
        font_scale,
        (255, 255, 255),
        1,
        cv2.LINE_AA,
    )
    return img


def find_image_id(filename: str, coco_images: List[dict]) -> Optional[int]:
    """Match a local filename to a COCO image entry (case-insensitive,
    stem fallback)."""
    path = Path(filename)
    name = path.name.lower()
    stem = path.stem.lower()
    # handles double extension e.g. .jpg.jpg
    bare_stem = Path(stem).stem.lower()
    for img in coco_images:
        coco_path = Path(img["file_name"])
        coco_name = coco_path.name.lower()
        coco_stem = coco_path.stem.lower()
        if coco_name == name or coco_stem == stem or coco_stem == bare_stem:
            return img["id"]
    return None


def main():
    args = parse_arguments()

    show_bbox = _bool(args.show_bbox)
    show_mask = _bool(args.show_mask)
    show_labels = _bool(args.show_labels)
    show_count = _bool(args.show_count)

    # Load COCO
    print(f"Loading annotations: {args.annotations}")
    with open(args.annotations) as f:
        coco_data = json.load(f)

    category_map: Dict[int, str] = {
        cat["id"]: cat["name"] for cat in coco_data.get("categories", [])
    }
    print(f"Categories found: {list(category_map.values())}")

    # Category filter
    allowed_cat_ids: Optional[set] = None
    if args.filter_categories.strip():
        wanted = {n.strip().lower() for n in args.filter_categories.split(",")}
        allowed_cat_ids = {
            cid for cid, name in category_map.items() if name.lower() in wanted
        }
        matched = {category_map[cid] for cid in allowed_cat_ids}
        print(f"Filtering to: {matched}")
        if not allowed_cat_ids:
            print("Warning: no categories matched — showing all.")
            allowed_cat_ids = None

    category_colors = (
        {cid: color_for_id(i) for i, cid in enumerate(sorted(category_map))}
        if args.color_mode == "per_category"
        else {}
    )

    # Group annotations by image_id, applying filters
    ann_by_image: Dict[int, List[dict]] = {}
    for ann in coco_data.get("annotations", []):
        if (
            allowed_cat_ids is not None
            and ann.get("category_id") not in allowed_cat_ids
        ):
            continue
        ann_by_image.setdefault(ann["image_id"], []).append(ann)

    output_mode = args.output_mode

    # Output dirs / paths
    frames_dir = Path(args.outdir) / "outputs_annotated"
    video_path = Path(args.outdir) / "output_video.mp4"
    if output_mode in ("frames", "both"):
        frames_dir.mkdir(parents=True, exist_ok=True)
    else:
        Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # Discover input images
    data_dir = Path("data_files")
    image_files = sorted(
        f
        for f in data_dir.glob("*")
        if f.suffix.lower() in {".jpg", ".jpeg", ".png", ".tiff", ".tif"}
    )

    if not image_files:
        print("Error: no image files found in data_files/")
        return

    print(f"\nProcessing {len(image_files)} image(s)...")
    matched_count = 0
    unmatched_count = 0
    video_frames: List[np.ndarray] = []

    for img_path in image_files:
        image = cv2.imread(str(img_path))
        if image is None:
            print(f"  Error: could not read {img_path.name}")
            continue

        image_id = find_image_id(img_path.name, coco_data.get("images", []))
        if image_id is None:
            print(
                f"{img_path.name}: no COCO match — output without annotations"
            )
            unmatched_count += 1
            continue
        else:
            anns = ann_by_image.get(image_id, [])
            print(f"  {img_path.name}: {len(anns)} annotation(s)")
            matched_count += 1

        annotated = draw_annotations(
            image,
            anns,
            category_map,
            show_bbox=show_bbox,
            show_mask=show_mask,
            show_labels=show_labels,
            mask_opacity=args.mask_opacity,
            bbox_thickness=args.bbox_thickness,
            color_mode=args.color_mode,
            font_scale=args.font_scale,
            category_colors=category_colors,
        )

        if show_count:
            draw_count_overlay(annotated, len(anns), args.font_scale)

        if output_mode in ("frames", "both"):
            out_path = (
               frames_dir / f"{Path(img_path.stem).stem}.{args.output_format}"
            )
            print(f"out_path {out_path}")
            cv2.imwrite(str(out_path), annotated)

        if output_mode in ("video", "both"):
            video_frames.append(annotated)

    # Compile video
    if output_mode in ("video", "both") and video_frames:
        h, w = video_frames[0].shape[:2]
        fourcc = cv2.VideoWriter_fourcc(*"mp4v")
        writer = cv2.VideoWriter(
            str(video_path), fourcc, args.video_fps, (w, h)
        )
        for frame in video_frames:
            writer.write(frame)
        writer.release()
        print(f"  Video: {video_path}")

    print(f"\n✓ Done: {matched_count} annotated, {unmatched_count} unmatched")
    if output_mode in ("frames", "both"):
        print(f"  Frames: {frames_dir}")
    if output_mode in ("video", "both"):
        print(f"  Video:  {video_path}")


if __name__ == "__main__":
    main()
