#!/usr/bin/env python3

import argparse
import json
import re
import sys
from pathlib import Path

from langchain_text_splitters import (
    CharacterTextSplitter,
    NLTKTextSplitter,
    RecursiveCharacterTextSplitter,
    SpacyTextSplitter,
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

# Reasons why a chunk cannot be traced back to the input, see diagnose_chunk().
CHUNK_TEXT_ALTERED = "chunk_text_altered"
START_INDEX_INVALID = "start_index_invalid"


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
        choices=("recursive_character", "character", "nltk", "spacy", "token"),
        default="recursive_character",
    )
    parser.add_argument("--sentence-language", default="english")
    parser.add_argument(
        "--spacy-pipeline",
        choices=("sentencizer", "en_core_web_sm"),
        default="sentencizer",
    )
    parser.add_argument(
        "--spacy-max-length",
        type=int,
        default=1_000_000,
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
        text = input_path.read_bytes().decode("utf-8")
    except UnicodeDecodeError as error:
        sys.exit(
            "The input dataset is not valid UTF-8 text "
            f"(invalid byte at offset {error.start}). "
            "Convert the dataset to UTF-8 before splitting it."
        )

    # Galaxy appends a newline to an uploaded file that has none, see
    # convert_newlines() in galaxy/lib/galaxy/datatypes/sniff.py. It is
    # deliberately kept here. An uploaded "abc" and an uploaded "abc\n" both
    # arrive as "abc\n", so the two cases cannot be told apart, and removing the
    # newline would cut a real character off every input that legitimately ends
    # in one, which is the common case. The chunks would then no longer add up
    # to the input, and the character based splitters would become lossy for no
    # gain. It would not even settle the NLTK case, because punkt discards all
    # trailing whitespace, so an input ending in a blank line still loses more
    # than the one appended character. What NLTK drops is reported instead, see
    # the warning at the end of main().
    return text


def get_tiktoken_options(args):
    allow_all_special = args.allowed_special == "all"

    return {
        "encoding_name": args.encoding_name,
        "model_name": args.model_name or None,
        "allowed_special": "all" if allow_all_special else set(),
        "disallowed_special": () if allow_all_special else "all",
    }


def special_token_error(error):
    """Return the error to raise for a ValueError coming from tiktoken.

    tiktoken raises when the input contains a special token string such as
    <|endoftext|> while the user asked for those strings to be rejected. Any
    other ValueError is handed back unchanged.
    """
    if "disallowed special token" not in str(error):
        return error

    return SystemExit(
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
            raise special_token_error(error)

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

    length_options = {
        **common_options,
        "length_function": length_function,
    }

    if args.splitter_type == "nltk":
        # strip_whitespace is handled on the finished chunks in main() instead of
        # here. separator="" together with strip_whitespace=False keeps the chunk text
        # identical to the matching slice of the input, so the reported start
        # indices stay usable. Note that the span based tokenizer drops whatever
        # follows the last sentence, see the warning in main().
        return NLTKTextSplitter(
            **length_options,
            language=args.sentence_language,
            separator="",
            use_span_tokenize=True,
            strip_whitespace=False,
        )

    if args.splitter_type == "spacy":
        # strip_whitespace is handled on the finished chunks in main() instead of
        # here. langchain strips every sentence before joining them, and
        # separator="" joins with nothing, so letting it strip would also delete
        # the space *between* two sentences and produce "mat.The".
        splitter = SpacyTextSplitter(
            **length_options,
            pipeline=args.spacy_pipeline,
            max_length=args.spacy_max_length,
            separator="",
            strip_whitespace=False,
        )
        # langchain only forwards max_length to the pipelines it loads with
        # spacy.load(). The sentencizer is built from English() instead and
        # silently keeps the spaCy default of one million characters, so it is
        # applied here for both. Guarded because _tokenizer is private API.
        tokenizer = getattr(splitter, "_tokenizer", None)

        if hasattr(tokenizer, "max_length"):
            tokenizer.max_length = args.spacy_max_length

        return splitter

    character_options = {
        **length_options,
        "strip_whitespace": args.strip_whitespace,
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


def diagnose_chunk(
    start_index,
    chunk_text,
    input_text,
    previous_start_index,
):
    """Return None when the chunk and its reported position are sound.

    Otherwise return the reason, so that the two very different causes can be
    reported separately: text the splitter altered, and a position the splitter
    got wrong.
    """
    if (
        start_index >= 0
        and (previous_start_index is None or start_index > previous_start_index)
        and input_text.startswith(chunk_text, start_index)
    ):
        return None

    # Only reached when something is already wrong, so scanning the whole input
    # once more is acceptable here.
    if chunk_text not in input_text:
        return CHUNK_TEXT_ALTERED

    return START_INDEX_INVALID


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

    # Checked before the tokenizer is built so that a job that cannot run does
    # not first pay for loading the tiktoken encoding.
    if args.splitter_type == "spacy" and len(input_text) > args.spacy_max_length:
        sys.exit(
            f"The input dataset holds {len(input_text)} characters, which is "
            "more than the configured spaCy maximum input length of "
            f"{args.spacy_max_length}. Raise 'Maximum input length' or split "
            "the dataset into smaller parts first. Mind the memory cost per "
            "input character stated in the help of that setting."
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

    try:
        splitter = build_splitter(
            args,
            input_text,
            length_function,
            tiktoken_options,
        )
    except OSError as error:
        if args.splitter_type == "spacy":
            sys.exit(
                f"The selected spaCy pipeline {args.spacy_pipeline!r} "
                "is not installed in the tool environment.\n"
                f"Original error: {error}"
            )
        raise
    except LookupError as error:
        # NLTKTextSplitter loads the punkt_tab data in its constructor, so this
        # has to be caught around build_splitter() and not around
        # create_documents().
        # The data comes from the nltk_data requirement, so this is reached only
        # on an environment that does not provide it, or provides it outside the
        # directories nltk searches:
        # -/usr/share/nltk_data
        # -/usr/local/share/nltk_data
        # -/usr/lib/nltk_data
        # -/usr/local/lib/nltk_data
        # -any directory listed in NLTK_DATA
        # KeyError and IndexError also derive from LookupError, so they are
        # excluded to keep an unrelated failure from being reported as missing
        # Punkt data.
        if args.splitter_type == "nltk" and not isinstance(
            error, (KeyError, IndexError)
        ):
            sys.exit(
                "The NLTK Punkt data required for the selected language is not "
                "installed. The Galaxy tool environment must provide the "
                "'punkt_tab' NLTK resource, for example below a directory "
                "listed in the NLTK_DATA environment variable. Use the spaCy "
                "sentence splitter if the data cannot be installed.\n"
                f"Original error: {error}"
            )
        raise

    try:
        documents = splitter.create_documents([input_text])
    except ValueError as error:
        raise special_token_error(error)

    chunks = []
    start_index_warnings = []
    altered_text_warnings = []
    previous_start_index = None
    empty_chunks = 0

    for document in documents:
        raw_chunk_text = document.page_content
        start_index = document.metadata.get("start_index")

        if args.strip_whitespace:
            # Only the outer whitespace of the finished chunk is removed, which
            # is what the option promises. The chunk therefore stays a verbatim
            # slice of the input and its start index stays valid, so it only has
            # to move by whatever came off the front.
            #
            # The splitter has already packed the chunks at this point, counting
            # the whitespace that is removed here, so a stripped chunk can come
            # out shorter than the requested size. That is the conservative
            # direction, and measuring the stripped text instead would mean
            # stripping before the packing, which is exactly what produces
            # "mat.The".
            leading = len(raw_chunk_text) - len(raw_chunk_text.lstrip())
            raw_chunk_text = raw_chunk_text.strip()

            if start_index is not None and start_index >= 0:
                start_index += leading

            # Nothing is left for a downstream step to work on. Reachable only
            # when the user asked for stripping, so no content is lost: without
            # it every chunk is kept, whitespace included.
            if not raw_chunk_text:
                empty_chunks += 1
                continue

        chunk_number = len(chunks) + 1

        # TODO: report the invalid start index upstream to langchain-text-splitters, since this
        # should not happen. Note that the invalid indices are not always negative: once one chunk
        # gets a negative index, the positions derived from it stay positive but point at the wrong
        # place, which is why diagnose_chunk() also checks the order and the content.
        # Potential cause: langchain text splitters might calculate the start index based on the length of the chunk in characters, even though the chunk size is measured in tokens.
        # reproduce with Character splitter recursive
        # input: test-data/langchain_docs_sample.txt
        # Text splitting separators and their order of use: default
        # Place separator at start of the following chunk start → ["First paragraph.", "\n\nSecond paragraph."]
        # Strip whitespace around chunks: true
        # Target chunk size should be determined by the number of tokens
        # Target chunk size: 100
        # Chunk overlap: 20
        problem = (
            None
            if start_index is None
            else diagnose_chunk(
                start_index,
                raw_chunk_text,
                input_text,
                previous_start_index,
            )
        )

        if problem is None:
            if start_index is not None:
                previous_start_index = start_index
        else:
            # The position is reported as null so that the JSON field keeps a
            # single type. The warnings below carry the detail.
            if problem == CHUNK_TEXT_ALTERED:
                altered_text_warnings.append(chunk_number)
            else:
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

    if altered_text_warnings:
        print(
            "WARNING: The text of the following chunk(s) does not occur in the "
            f"input: {', '.join(str(number) for number in altered_text_warnings)}. "
            "The splitter returned text it had modified, so no position in the "
            "input describes it and the start index is reported as null. The "
            "known cause is splitting between tokens: a cut can fall inside a "
            "character that is encoded in several bytes, which then becomes the "
            "Unicode replacement character. Use one of the character based "
            "splitters for text that is not plain ASCII.",
            flush=True,
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

    if empty_chunks:
        print(
            f"WARNING: {empty_chunks} chunk(s) held nothing but whitespace and "
            "were left out of the outputs, because stripping the whitespace "
            "around the chunks left them empty.",
            flush=True,
        )

    # The NLTK splitter builds the chunks from the sentence spans reported by
    # punkt, which end at the last sentence. Punkt trims trailing whitespace, so
    # whatever follows is dropped and the chunks no longer add up to the input.
    # The first span always starts at offset 0, so nothing is lost in front.
    # The spaCy splitter keeps everything.
    if args.splitter_type == "nltk" and chunks and not args.strip_whitespace:
        last_chunk = chunks[-1]
        last_text = last_chunk["text"]
        # The validated start index is preferred over searching for the text:
        # rfind() returns the *last* occurrence, so it would report no loss
        # whenever the dropped text happens to repeat the final chunk.
        last_start = last_chunk["start_index"]

        if last_start is None:
            last_start = input_text.rfind(last_text)

        dropped = (
            len(input_text) - (last_start + len(last_text)) if last_start >= 0 else 0
        )

        if dropped > 0:
            print(
                f"WARNING: The NLTK sentence splitter dropped the {dropped} "
                "character(s) after the last sentence. Use the spaCy sentence "
                "splitter to keep the complete input.",
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
        # Reported from the argument rather than per splitter, because the
        # stripping is done above for whichever splitter produced the chunks.
        # A splitter that does not offer the option never receives it, so the
        # value stays false for it without having to be hardcoded here.
        "strip_whitespace": args.strip_whitespace,
    }

    if args.splitter_type in ("character", "recursive_character"):
        metadata["keep_separator"] = args.keep_separator

    elif args.splitter_type == "nltk":
        metadata["sentence_language"] = args.sentence_language

    elif args.splitter_type == "spacy":
        metadata["spacy_pipeline"] = args.spacy_pipeline

    args.output_json.write_text(
        json.dumps(
            {**metadata, "chunks": chunks},
            indent=2,
        )
        + "\n",
        encoding="utf-8",
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
