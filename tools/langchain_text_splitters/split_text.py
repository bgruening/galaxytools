#!/usr/bin/env python3

import argparse
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

MAX_UNICODE_CODE_POINT = 0x10FFFF
SURROGATE_RANGE = range(0xD800, 0xE000)

# The chunk text is stored in a single TSV cell, so every character that would
# otherwise start a new row or a new column has to be escaped. The backslash is
# escaped first so that the transformation stays reversible.
TSV_ESCAPES = (
    ("\\", r"\\"),
    ("\t", r"\t"),
    ("\r", r"\r"),
    ("\n", r"\n"),
)

LENGTH_KEY = "length"


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
    tokenizer_group = parser.add_mutually_exclusive_group()
    tokenizer_group.add_argument("--encoding-name", default="gpt2")
    tokenizer_group.add_argument("--model-name", default="")
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


def read_input_text(input_path):
    # The bytes are decoded explicitly instead of using Path.read_text(), which
    # would apply universal newline translation and silently rewrite \r\n and \r
    # line endings as \n. That would change both the chunk content and the
    # reported start indices with respect to the input dataset.
    try:
        return input_path.read_bytes().decode("utf-8")
    except UnicodeDecodeError as error:
        sys.exit(
            "The input dataset is not valid UTF-8 text "
            f"(invalid byte at offset {error.start}). "
            "Convert the dataset to UTF-8 before splitting it."
        )


def get_tiktoken_options(args):
    allow_all_special = args.allowed_special == "all"

    return {
        "encoding_name": args.encoding_name,
        "model_name": args.model_name or None,
        "allowed_special": "all" if allow_all_special else set(),
        "disallowed_special": () if allow_all_special else "all",
    }


def fail_on_disallowed_special_token(error):
    # Raised by tiktoken when the input contains a special token string such as
    # <|endoftext|> while the user asked for those strings to be rejected.
    if "disallowed special token" not in str(error):
        raise error

    sys.exit(
        "The input contains a special token string that is rejected by the "
        "selected tokenizer. Set the special token handling to allow all "
        "special token strings, or remove the special token from the input.\n"
        f"Original error: {error}"
    )


def build_tiktoken_counter(
    encoding_name,
    model_name=None,
    allowed_special=None,
    disallowed_special="all",
):
    import tiktoken

    if allowed_special is None:
        allowed_special = set()

    # Note: below needs either internet access to download the encodings from
    # the OpenAI servers, or the TIKTOKEN_CACHE_DIR environment variable
    # pointing to a directory that already holds the encoding files. The
    # tiktoken conda package does not ship them.
    if model_name:
        encoding = tiktoken.encoding_for_model(model_name)
    else:
        encoding = tiktoken.get_encoding(encoding_name)

    def count_tokens(text):
        try:
            tokens = encoding.encode(
                text,
                allowed_special=allowed_special,
                disallowed_special=disallowed_special,
            )
        except ValueError as error:
            fail_on_disallowed_special_token(error)
            raise

        return len(tokens)

    return count_tokens


def decode_custom_separator(value):
    simple_escapes = {
        r"\n": "\n",
        r"\r": "\r",
        r"\t": "\t",
        r"\\": "\\",
    }

    decoded = []
    position = 0

    while position < len(value):
        character = value[position]

        if character != "\\":
            decoded.append(character)
            position += 1
            continue

        match = CUSTOM_SEPARATOR_ESCAPE_PATTERN.match(value, position)

        if match is None:
            sys.exit(
                "The custom separator contains an unsupported escape sequence "
                f"at position {position}: {value[position:position + 10]!r}. "
                r"Supported escapes are \\, \n, \r, \t, "
                r"\uXXXX (4 hex digits) and \UXXXXXXXX (8 hex digits)."
            )

        escape = match.group(0)

        if escape in simple_escapes:
            decoded.append(simple_escapes[escape])
        else:
            code_point = int(escape[2:], 16)

            if code_point > MAX_UNICODE_CODE_POINT or code_point in SURROGATE_RANGE:
                sys.exit(
                    f"The custom separator escape {escape} is not a valid "
                    "Unicode character."
                )

            decoded.append(chr(code_point))

        position = match.end()

    return "".join(decoded)


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
        separators = resolve_separator_specs(args.separator_specs)
        # Append an empty string as the final fallback separator to ensure that the text can always be split,
        # even if none of the tried separators before were able to split the text without exceeding the chunk size.
        character_options["separators"] = [*separators, ""]

    return RecursiveCharacterTextSplitter(**character_options)


def write_text_output(
    output_path,
    metadata,
    chunks,
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
                    f"{chunk[LENGTH_KEY]} {length_label}, "
                    f"start={chunk['start_index']} ---"
                ),
                chunk["text"],
                "",
            ]
        )
    output_path.write_text("\n".join(lines), encoding="utf-8")


def write_chunk_files(chunks_dir, chunks):
    chunks_dir.mkdir(parents=True, exist_ok=True)

    for chunk in chunks:
        chunk_path = chunks_dir / f"chunk_{chunk['index']:04d}.txt"
        chunk_path.write_text(
            chunk["text"],
            encoding="utf-8",
        )


def escape_tsv_text(text):
    for raw, escaped in TSV_ESCAPES:
        text = text.replace(raw, escaped)

    return text


def write_tsv_output(output_path, chunks):
    # The rows are written without the csv module on purpose. Its writer would
    # additionally apply CSV quoting to any chunk containing a double quote,
    # which the documented escaping above cannot undo. escape_tsv_text() already
    # removes every character that could break the column or row structure.
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        for chunk in chunks:
            fields = [
                str(chunk["index"]),
                escape_tsv_text(chunk["text"]),
                str(chunk[LENGTH_KEY]),
            ]
            handle.write("\t".join(fields) + "\n")


def main():
    args = parse_args()

    if args.chunk_overlap >= args.chunk_size:
        sys.exit("Chunk overlap must be smaller than chunk size.")

    input_text = read_input_text(args.input)

    if not input_text.strip():
        sys.exit(
            "The input dataset is empty or contains only whitespace. "
            "There is nothing to split."
        )

    length_mode = "token" if args.splitter_type == "token" else args.length_mode

    tiktoken_options = get_tiktoken_options(args)

    if length_mode == "token":
        length_function = build_tiktoken_counter(**tiktoken_options)
        length_label = "tokens"
        tokenizer_label = args.model_name or args.encoding_name
        length_function_label = f"tiktoken:{tokenizer_label}"
    else:
        length_function = len
        length_label = "characters"
        length_function_label = "characters"

    splitter = build_splitter(
        args,
        input_text,
        length_function,
        tiktoken_options,
    )

    try:
        documents = splitter.create_documents([input_text])
    except ValueError as error:
        fail_on_disallowed_special_token(error)
        raise

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
            # The invalid index is reported as null so that the JSON keeps a
            # single type for the field. The warning below carries the detail.
            start_index_warnings.append((chunk_number, start_index))
            start_index = None

        chunks.append(
            {
                "index": chunk_number,
                LENGTH_KEY: length_function(raw_chunk_text),
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
            f"invalid start index returned for chunk(s): {affected_chunks}. "
            "The start index of these chunks is reported as null.",
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

    if args.splitter_type == "token":
        # The token splitter has no separators and never strips whitespace, so
        # the value is reported as false regardless of what was requested.
        metadata["strip_whitespace"] = False
    else:
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
        length_label,
    )
    write_tsv_output(
        args.output_tsv,
        chunks,
    )
    write_chunk_files(
        args.chunks_dir,
        chunks,
    )


if __name__ == "__main__":
    main()
