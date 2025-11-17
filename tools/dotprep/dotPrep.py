#! /usr/bin/env python

# Author: Maria Nattestad
# Email: maria.nattestad@gmail.com

import argparse
import gzip
import operator
import re
import time


import numpy as np


def run(args):
    filename = args.delta
    unique_length = args.unique_length
    output_filename = args.out
    keep_small_uniques = True
    max_overview_alignments = args.overview

    header_lines_by_query, lines_by_query = get_query_ref_combinations(filename)

    unique_alignments = calculate_uniqueness(
        header_lines_by_query, lines_by_query, unique_length, keep_small_uniques
    )

    reference_lengths, fields_by_query = write_filtered_delta_file(
        filename,
        output_filename,
        unique_alignments,
        unique_length,
        header_lines_by_query,
    )

    index_for_dot(reference_lengths, fields_by_query, output_filename, max_overview_alignments)


def scrub(string):
    return (
        string.replace(",", "_")
        .replace("!", "_")
        .replace("~", "_")
        .replace("#", "_")
    )


def get_query_ref_combinations(filename):
    print("header from delta file:")
    try:
        f = gzip.open(filename, "rt")
        print(f.readline().strip())
    except OSError:
        f = open(filename, "r")
        print(f.readline().strip())

    print(f.readline().strip())
    linecounter = 0
    current_query_name = ""
    current_header = ""
    lines_by_query = {}
    header_lines_by_query = {}
    before = time.time()

    for line in f:
        if line[0] == ">":
            linecounter += 1
            current_header = line.strip()
            current_query_name = scrub(current_header.split()[1])

            if header_lines_by_query.get(current_query_name) is None:
                lines_by_query[current_query_name] = []
                header_lines_by_query[current_query_name] = []
        else:
            fields = line.strip().split()
            if len(fields) > 4:
                query_min = min(int(fields[2]), int(fields[3]))
                query_max = max(int(fields[2]), int(fields[3]))
                lines_by_query[current_query_name].append((query_min, query_max))
                header_lines_by_query[current_query_name].append(current_header)

    f.close()
    print(
        "First read through the file: %d seconds for %d query-reference combinations"
        % (time.time() - before, linecounter)
    )
    return header_lines_by_query, lines_by_query


def calculate_uniqueness(header_lines_by_query, lines_by_query, unique_length, keep_small_uniques):
    before = time.time()
    unique_alignments = {}
    num_queries = len(lines_by_query)
    print("Filtering alignments of %d queries" % num_queries)

    num_query_step_to_report = max(num_queries / 100, 1)
    query_counter = 0

    for query in lines_by_query:
        unique_alignments[query] = summarize_planesweep(
            lines_by_query[query],
            unique_length_required=unique_length,
            keep_small_uniques=keep_small_uniques,
        )
        query_counter += 1
        if query_counter % num_query_step_to_report == 0:
            print("Progress: %d%%" % (query_counter * 100 / num_queries))

    print("Progress: 100%")
    print(
        "Deciding which alignments to keep: %d seconds for %d queries"
        % (time.time() - before, num_queries)
    )
    return unique_alignments


def summarize_planesweep(lines, unique_length_required, keep_small_uniques=False):
    unique_alignments = []
    if len(lines) == 0:
        return []
    if len(lines) == 1:
        if keep_small_uniques or abs(lines[0][1] - lines[0][0]) >= unique_length_required:
            return [0]
        return []

    starts_and_stops = []
    for query_min, query_max in lines:
        starts_and_stops.append((query_min, "start"))
        starts_and_stops.append((query_max, "stop"))

    sorted_starts_and_stops = sorted(starts_and_stops, key=operator.itemgetter(0))
    current_coverage = 0
    last_position = -1
    sorted_unique_intervals_left = []
    sorted_unique_intervals_right = []

    for pos, change in sorted_starts_and_stops:
        if current_coverage == 1:
            sorted_unique_intervals_left.append(last_position)
            sorted_unique_intervals_right.append(pos)
        if change == "start":
            current_coverage += 1
        else:
            current_coverage -= 1
        last_position = pos

    linecounter = 0
    for query_min, query_max in lines:
        i = binary_search(query_min, sorted_unique_intervals_left, 0, len(sorted_unique_intervals_left))
        exact_match = False
        if i < len(sorted_unique_intervals_left):
            if (
                sorted_unique_intervals_left[i] == query_min
                and sorted_unique_intervals_right[i] == query_max
            ):
                exact_match = True
        sum_uniq = 0
        while (
            i < len(sorted_unique_intervals_left)
            and sorted_unique_intervals_left[i] >= query_min
            and sorted_unique_intervals_right[i] <= query_max
        ):
            sum_uniq += sorted_unique_intervals_right[i] - sorted_unique_intervals_left[i]
            i += 1
        if sum_uniq >= unique_length_required:
            unique_alignments.append(linecounter)
        elif keep_small_uniques and exact_match:
            unique_alignments.append(linecounter)
        linecounter += 1

    return unique_alignments


def binary_search(query, numbers, left, right):
    if left >= right:
        return right
    mid = int((right + left) / 2)
    if query == numbers[mid]:
        return mid
    if query < numbers[mid]:
        return binary_search(query, numbers, left, mid)
    return binary_search(query, numbers, mid + 1, right)


def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r"(\d+)", string_)]


def write_filtered_delta_file(filename, output_filename, unique_alignments, unique_length, header_lines_by_query):
    before = time.time()
    f_out_delta = gzip.open(
        f"{output_filename}.uniqueAnchorFiltered_l{unique_length}.delta.gz",
        "wt",
    )
    try:
        f = gzip.open(filename, "rt")
        header1 = f.readline()
    except OSError:
        f = open(filename, "r")
        header1 = f.readline()

    f_out_delta.write(header1)
    f_out_delta.write(f.readline())
    linecounter = 0

    list_of_unique_alignments = []
    alignment_counter = {}
    keep_printing = False
    ref_sequences = set()
    query_sequences = set()
    reference_lengths = []
    query_lengths = {}
    fields_by_query = {}

    for line in f:
        linecounter += 1
        if line[0] == ">":
            fields = line.strip().split()
            query = scrub(fields[1])
            list_of_unique_alignments = unique_alignments[query]

            header_needed = any(
                line.strip() == header_lines_by_query[query][index]
                for index in list_of_unique_alignments
            )
            if header_needed:
                f_out_delta.write(line)

            alignment_counter[query] = alignment_counter.get(query, 0)
            current_reference_name = scrub(fields[0][1:])
            current_query_name = scrub(fields[1])
            current_reference_size = int(fields[2])
            current_query_size = int(fields[3])

            if current_reference_name not in ref_sequences:
                reference_lengths.append((current_reference_name, current_reference_size))
                ref_sequences.add(current_reference_name)
            if current_query_name not in query_sequences:
                query_lengths[current_query_name] = current_query_size
                query_sequences.add(current_query_name)

        else:
            fields = line.strip().split()
            if len(fields) > 4:
                ref_start = int(fields[0])
                ref_end = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                csv_tag = "repetitive"
                if alignment_counter[query] in list_of_unique_alignments:
                    f_out_delta.write(line)
                    csv_tag = "unique"
                    keep_printing = True
                else:
                    keep_printing = False
                record = [
                    ref_start,
                    ref_end,
                    query_start,
                    query_end,
                    current_reference_size,
                    current_query_size,
                    current_reference_name,
                    current_query_name,
                    csv_tag,
                ]
                if fields_by_query.get(current_query_name) is None:
                    fields_by_query[current_query_name] = []
                fields_by_query[current_query_name].append(record)
                alignment_counter[query] += 1
            elif keep_printing:
                f_out_delta.write(line)

    f.close()
    f_out_delta.close()
    print(
        "Writing filtered delta file and capturing information for coords file: "
        "%d seconds for %d total lines in file"
        % (time.time() - before, linecounter)
    )
    return reference_lengths, fields_by_query


def index_for_dot(reference_lengths, fields_by_query, output_prefix, max_overview_alignments):
    reference_lengths.sort(key=lambda x: natural_key(x[0]))
    cumulative_sum = 0
    ref_chrom_offsets = {}
    queries_by_reference = {}

    for ref_name, ref_length in reference_lengths:
        ref_chrom_offsets[ref_name] = cumulative_sum
        cumulative_sum += ref_length
        queries_by_reference[ref_name] = set()

    flip_by_query = {}
    unique_references_by_query = {}
    all_references_by_query = {}
    relative_ref_position_by_query = []
    ordered_tags = ["unique", "repetitive"]

    f_out_coords = open(f"{output_prefix}.coords", "w", encoding="utf-8")
    f_out_coords.write("ref_start,ref_end,query_start,query_end,ref\n")

    query_byte_positions = {}
    query_lengths = {}
    all_alignments = []
    last_query = ""

    for query_name, lines in fields_by_query.items():
        sum_forward = 0
        sum_reverse = 0
        ref_position_scores = []
        unique_references_by_query[query_name] = set()
        all_references_by_query[query_name] = set()

        for fields in lines:
            tag = fields[8]
            query_lengths[query_name] = int(fields[5])
            ref_name = fields[6]
            all_references_by_query[query_name].add(ref_name)
            if tag == "unique":
                query_start = int(fields[2])
                query_stop = int(fields[3])
                ref_start = int(fields[0])
                ref_stop = int(fields[1])
                alignment_length = abs(query_stop - query_start)
                unique_references_by_query[query_name].add(ref_name)
                queries_by_reference[ref_name].add(query_name)
                ref_position_scores.append(
                    ref_chrom_offsets[ref_name] + (ref_start + ref_stop) / 2
                )
                if query_stop < query_start:
                    sum_reverse += alignment_length
                else:
                    sum_forward += alignment_length

        flip = sum_reverse > sum_forward
        flip_by_query[query_name] = "-" if flip else "+"

        for tag in ordered_tags:
            query_byte_positions[(last_query, "end")] = f_out_coords.tell()
            query_byte_positions[(query_name, tag)] = f_out_coords.tell()
            f_out_coords.write(f"!{query_name}!{tag}\n")
            for fields in lines:
                if fields[8] == tag:
                    if flip:
                        fields[2] = int(fields[5]) - int(fields[2])
                        fields[3] = int(fields[5]) - int(fields[3])
                    output_fields = [fields[0], fields[1], fields[2], fields[3], fields[6]]
                    f_out_coords.write(",".join(str(i) for i in output_fields) + "\n")
                    alignment_length = abs(int(fields[3]) - int(fields[2]))
                    all_alignments.append((fields, alignment_length))

        rel_pos = np.median(ref_position_scores) if ref_position_scores else 0
        relative_ref_position_by_query.append((query_name, rel_pos))
        last_query = query_name

    query_byte_positions[(last_query, "end")] = f_out_coords.tell()
    relative_ref_position_by_query.sort(key=lambda x: x[1])

    f_out_index = open(f"{output_prefix}.coords.idx", "w", encoding="utf-8")
    f_out_index.write("#ref\n")
    f_out_index.write("ref,ref_length,matching_queries\n")
    for ref_name, ref_length in reference_lengths:
        f_out_index.write(
            f"{ref_name},{ref_length},{'~'.join(queries_by_reference[ref_name])}\n"
        )

    f_out_index.write("#query\n")
    f_out_index.write(
        "query,query_length,orientation,bytePosition_unique,"
        "bytePosition_repetitive,bytePosition_end,unique_matching_refs,matching_refs\n"
    )

    for query_name, _ in relative_ref_position_by_query:
        f_out_index.write(
            f"{query_name},{query_lengths[query_name]},{flip_by_query[query_name]},"
            f"{query_byte_positions[(query_name, 'unique')]},"
            f"{query_byte_positions[(query_name, 'repetitive')] - query_byte_positions[(query_name, 'unique')]},"
            f"{query_byte_positions[(query_name, 'end')] - query_byte_positions[(query_name, 'repetitive')]},"
            f"{'~'.join(unique_references_by_query[query_name])},"
            f"{'~'.join(all_references_by_query[query_name])}\n"
        )

    f_out_index.write("#overview\n")
    f_out_index.write("ref_start,ref_end,query_start,query_end,ref,query,tag\n")

    num_overview_alignments = min(max_overview_alignments, len(all_alignments))
    if num_overview_alignments < len(all_alignments):
        print(
            f"Included the longest {max_overview_alignments} alignments in the index "
            f"under #overview (change this with the --overview parameter), "
            f"out of a total of {len(all_alignments)} alignments."
        )

    all_alignments.sort(key=lambda x: -x[1])
    for fields, _ in all_alignments[:num_overview_alignments]:
        f_out_index.write(",".join(str(i) for i in fields[:9]) + "\n")

    f_out_coords.close()
    f_out_index.close()


def main():
    parser = argparse.ArgumentParser(
        description="Take a delta file, apply Assemblytics unique anchor filtering, "
        "and prepare coordinates input files for Dot"
    )
    parser.add_argument("--delta", help="delta file", dest="delta", type=str, required=True)
    parser.add_argument("--out", help="output file", dest="out", type=str, default="output")
    parser.add_argument(
        "--unique-length",
        help="Minimum unique query sequence length required (default: 10000)",
        dest="unique_length",
        type=int,
        default=10000,
    )
    parser.add_argument(
        "--overview",
        help="Number of alignments to include in coords.idx overview (default: 1000)",
        dest="overview",
        type=int,
        default=1000,
    )
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
