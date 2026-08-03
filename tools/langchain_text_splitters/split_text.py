#!/usr/bin/env python3

import argparse
import csv
import json
import re
import sys
from pathlib import Path

from langchain_text_splitters import (
    CharacterTextSplitter,
    RecursiveCharacterTextSplitter,
    TokenTextSplitter,
)

SEPARATOR_MAP = {
    "paragraph_break": "\n\n",
    "line_break": "\n",
    "tab": "\t",
    "space": " ",
    "ascii_full_stop": ".",
    "ascii_comma": ",",
    "ascii_semicolon": ";",
    "zero_width_space": "\u200b",
    "fullwidth_comma": "\uff0c",
    "ideographic_comma": "\u3001",
    "fullwidth_full_stop": "\uff0e",
    "ideographic_full_stop": "\u3002",
}

KEEP_SEPARATOR_VALUES = {
    "start": "start",
    "end": "end",
    "false": False,
}

CUSTOM_SEPARATOR_ESCAPE_PATTERN = re.compile(
    r"\\(?:[\\nrt]|u[0-9a-fA-F]{4}|U[0-9a-fA-F]{8})"
)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Split text with langchain-text-splitters."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-text", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--chunks-dir", type=Path, required=True)
    parser.add_argument(
        "--splitter-type",
        choices=("recursive_character", "character", "token"),
        default="recursive_character",
    )
    parser.add_argument("--chunk-size", type=int, required=True)
    parser.add_argument("--chunk-overlap", type=int, default=0)
    parser.add_argument(
        "--length-mode",
        choices=("characters", "token"),
        default="characters",
    )
    parser.add_argument("--encoding-name", default="gpt2")
    parser.add_argument("--model-name", default="")
    parser.add_argument(
        "--allowed-special",
        choices=("none", "all"),
        default="none",
    )
    parser.add_argument(
        "--keep-separator",
        choices=("start", "end", "false"),
        default="start",
    )
    parser.add_argument(
        "--separator-name",
        choices=SEPARATOR_MAP,
    )
    parser.add_argument("--separator")
    parser.add_argument("--separator-specs", nargs="+", default=None)
    parser.add_argument("--strip-whitespace", action="store_true")
    return parser.parse_args()


def get_tiktoken_options(args):
    allow_all_special = args.allowed_special == "all"

    return {
        "encoding_name": args.encoding_name,
        "model_name": args.model_name or None,
        "allowed_special": "all" if allow_all_special else set(),
        "disallowed_special": () if allow_all_special else "all",
    }


def build_tiktoken_counter(
    encoding_name,
    model_name=None,
    allowed_special=None,
    disallowed_special="all",
):
    import tiktoken

    if allowed_special is None:
        allowed_special = set()
    
    # Note: below needs either the TIKTOKEN_CACHE_DIR environment variable set
    # respectively the directory populated with `python3 get_encodings.py` in the folder dev_utils.
    # or internet access to download the encodings from the OpenAI servers.
    if model_name:
        encoding = tiktoken.encoding_for_model(model_name)
    else:
        encoding = tiktoken.get_encoding(encoding_name)

    def count_tokens(text):
        return len(
            encoding.encode(
                text,
                allowed_special=allowed_special,
                disallowed_special=disallowed_special,
            )
        )

    return count_tokens

def decode_custom_separator(value):
    """Decode supported escape notation without altering literal Unicode."""
    simple_escapes = {
        r"\n": "\n",
        r"\r": "\r",
        r"\t": "\t",
        r"\\": "\\",
    }

    def replace_escape(match):
        escape = match.group(0)

        if escape in simple_escapes:
            return simple_escapes[escape]

        return chr(int(escape[2:], 16))

    return CUSTOM_SEPARATOR_ESCAPE_PATTERN.sub(replace_escape, value)


def resolve_separator(
    separator_name,
    custom_separator,
):
    if separator_name is not None:
        return SEPARATOR_MAP[separator_name]

    if custom_separator is not None:
        return decode_custom_separator(custom_separator)

    raise ValueError("No separator was provided.")


def resolve_separator_specs(specs):
    separators = []

    for spec in specs:
        kind, value = spec.split(":", 1)

        if kind == "predefined":
            separators.append(SEPARATOR_MAP[value])
        else:
            separators.append(decode_custom_separator(value))

    return separators


def build_splitter(args, input_text, length_function, tiktoken_options):
    common_options = {
        "chunk_size": args.chunk_size,
        "chunk_overlap": args.chunk_overlap,
        "add_start_index": True,
    }

    if args.splitter_type == "token":
        return TokenTextSplitter(
            **common_options,
            **tiktoken_options,
        )

    character_options = {
        **common_options,
        "strip_whitespace": args.strip_whitespace,
        "length_function": length_function,
        "keep_separator": KEEP_SEPARATOR_VALUES[args.keep_separator],
    }

    if args.splitter_type == "character":
        separator = resolve_separator(
            separator_name=args.separator_name,
            custom_separator=args.separator,
        )

        if separator not in input_text:
            sys.exit(
                "The selected separator "
                f"{separator!r} was not found in the input text. "
                "Select a separator that occurs in the input or use the "
                "recursive character splitter."
         )

        return CharacterTextSplitter(
            **character_options,
            separator=separator,
        )

    if args.separator_specs is not None:
        separators = (
            resolve_separator_specs(args.separator_specs)
            if args.separator_specs is not None
            else None
        )
        # Append an empty string as the final fallback separator to ensure that the text can always be split,
        # even if none of the tried separators before were able to split the text without exceeding the chunk size.
        character_options["separators"] = [*separators, ""]

    return RecursiveCharacterTextSplitter(**character_options)


def write_text_output(
    output_path,
    metadata,
    chunks,
    count_key,
    length_label,
):
    lines = [
        f"Splitter: {metadata['splitter_type']}",
        f"Length function: {metadata['length_function']}",
        f"Input characters: {metadata['input_characters']}",
        f"Input length: {metadata['input_length']} {length_label}",
        f"Chunk size: {metadata['chunk_size']}",
        f"Chunk overlap: {metadata['chunk_overlap']}",
        f"Number of chunks: {metadata['number_of_chunks']}",
        "",
    ]

    for chunk in chunks:
        lines.extend(
            [
                (
                    f"--- Chunk {chunk['index']}: "
                    f"{chunk[count_key]} {length_label}, "
                    f"start={chunk['start_index']} ---"
                ),
                chunk["text"],
                "",
            ]
        )
    output_path.write_text("\n".join(lines))


def write_chunk_files(chunks_dir, chunks):
    chunks_dir.mkdir(parents=True, exist_ok=True)

    for chunk in chunks:
        chunk_path = chunks_dir / f"chunk_{chunk['index']:04d}.txt"
        chunk_path.write_text(
            chunk["text"],
        )


def write_tsv_output(output_path, chunks, count_key):
    with output_path.open("w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        for chunk in chunks:
            table_text = chunk["text"].replace("\t", r"\t").replace("\n", r"\n")
            writer.writerow(
                [
                    chunk["index"],
                    table_text,
                    chunk[count_key],
                ]
            )


def main():
    args = parse_args()

    if args.chunk_overlap >= args.chunk_size:
        sys.exit("Chunk overlap must be smaller than chunk size.")

    input_text = args.input.read_text()

    length_mode = "token" if args.splitter_type == "token" else args.length_mode

    tiktoken_options = get_tiktoken_options(args)

    if length_mode == "token":
        length_function = build_tiktoken_counter(**tiktoken_options)
        count_key = "token_count"
        length_label = "tokens"
        tokenizer_label = args.model_name or args.encoding_name
        length_function_label = f"tiktoken:{tokenizer_label}"
    else:
        length_function = len
        count_key = "character_count"
        length_label = "characters"
        length_function_label = "characters"

    splitter = build_splitter(
        args,
        input_text,
        length_function,
        tiktoken_options,
    )
    documents = splitter.create_documents([input_text])

    chunks = []
    start_index_warnings = []

    for document in documents:
        raw_chunk_text = document.page_content
        start_index = document.metadata.get("start_index")
        chunk_number = len(chunks) + 1

        if start_index is not None and start_index < 0:
            # TODO: report upstream to langchain-text-splitters as a bug, since this should not happen.
            # Potential cause: langchain text splitters might calculate the start index based on the length of the chunk in characters, even though the chunk size is measured in tokens.
            # reproduce with Character splitter recursive
            # input: test-data/langchain_docs_sample.txt
            # Text splitting separators and their order of use: default
            # Place separator at start of the following chunk start → ["First paragraph.", "\n\nSecond paragraph."]
            # Strip whitespace around chunks: true
            # Target chunk size should be determined by the number of tokens
            # Target chunk size: 100
            # Chunk overlap: 20
            start_index_warnings.append((chunk_number, start_index))
            start_index = (
                "Potential upstream langchain-text-splitters bug: "
                f"received invalid start index {start_index!r}"
            )

        chunks.append(
            {
                "index": chunk_number,
                count_key: length_function(raw_chunk_text),
                "start_index": start_index,
                "text": raw_chunk_text,
            }
        )

    if start_index_warnings:
        affected_chunks = ", ".join(
            f"{chunk_number} (start_index: {start_index!r})"
            for chunk_number, start_index in start_index_warnings
        )
        print(
            "WARNING: Potential upstream langchain-text-splitters bug: "
            f"invalid start index returned for chunk(s): {affected_chunks}",
            flush=True,
        )

    metadata = {
        "input_characters": len(input_text),
        "input_length": length_function(input_text),
        "length_unit": length_label,
        "length_function": length_function_label,
        "splitter_type": args.splitter_type,
        "chunk_size": args.chunk_size,
        "chunk_overlap": args.chunk_overlap,
        "number_of_chunks": len(chunks),
    }

    if args.splitter_type != "token":
        metadata.update(
            {
                "keep_separator": args.keep_separator,
                "strip_whitespace": args.strip_whitespace,
            }
        )

    args.output_json.write_text(
        json.dumps(
            {**metadata, "chunks": chunks},
            indent=2,
        )
        + "\n",
    )

    write_text_output(
        args.output_text,
        metadata,
        chunks,
        count_key,
        length_label,
    )
    write_tsv_output(
        args.output_tsv,
        chunks,
        count_key,
    )
    write_chunk_files(
        args.chunks_dir,
        chunks,
    )


if __name__ == "__main__":
    main()
