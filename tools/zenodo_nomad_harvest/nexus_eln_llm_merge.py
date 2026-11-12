#!/usr/bin/env python
"""Check an LLM-written ELN and merge the values that pass into the rule ELN.

The answer of an ``llm_hub`` job holds two YAML blocks: the filled ELN and the
evidence (a quote and a file name per filled field). Only the fields of
``eln_rules.yaml`` (``fill`` when empty in the rule ELN, and ``replace``) are taken.
A value is kept when its quote is found in a source file (the same texts that
``nexus_llm_input`` wrote) and it passes the field's checks; its form may be cleaned.
Every kept, dropped or cleaned value is listed in the report. Without a usable answer
the output is the rule ELN.
"""

import argparse
import copy
import hashlib
import json
import logging
import re
import sys
import unicodedata
from pathlib import Path

import yaml

from nexus_llm_input import (
    BASE_ELN_NAME, INSTRUCTIONS_NAME, TOOL_DIR, definition_of, example_values, has_path,
    load_instructions, set_path, source_files,
)

logger = logging.getLogger("nexus_eln_llm_merge")

NAME_FIELDS = ("Sample.name", "Sample.Sample_component.name")
FORMULA_FIELD = "Sample.Sample_component.chemical_formula"
ELEMENTS = set("""
H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br
Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho
Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es
Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
""".split())
GREEK = {
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon", "ζ": "zeta",
    "η": "eta", "θ": "theta", "κ": "kappa", "λ": "lambda", "μ": "mu", "ν": "nu", "ξ": "xi",
    "π": "pi", "ρ": "rho", "σ": "sigma", "τ": "tau", "φ": "phi", "χ": "chi", "ψ": "psi",
    "ω": "omega",
}
DASHES = "‐‑‒–—−"  # hyphens, en and em dash, minus sign


class Dropped(Exception):
    """A value that fails a check; the message is the reason in the report."""


# A ":KEY:" line of a Nanonis .sxm header, e.g. ":NANONIS_VERSION:"; a key alone is no evidence.
HEADER_KEY = re.compile(r"(?<!\w):[A-Za-z][\w>/ -]*:(?!\w)")


def norm(text):
    """Text as quotes are compared: NFKC form, dashes as `-`, every run of whitespace as one space.

    Some models write typographic dashes and spaces; NFKC: https://unicode.org/reports/tr15/
    """
    text = unicodedata.normalize("NFKC", str(text))
    text = re.sub(f"[{DASHES}]", "-", text)
    return re.sub(r"\s+", " ", text).strip()


def squash(text):
    """For "value is in the text": no spaces, dashes or case; Greek letters spelled out."""
    text = str(text).lower()
    for letter, word in GREEK.items():
        text = text.replace(letter, word)
    return re.sub(rf"[\s\-{DASHES}]", "", text)


def value_at(eln, dotted):
    node = eln
    for key in dotted.split("."):
        if not isinstance(node, dict) or key not in node:
            return None
        node = node[key]
    return node


def parse_answer(text):
    """(ELN, evidence) from the answer's YAML blocks; the ELN is None if there is none."""
    eln, evidence = None, {}
    for block in re.findall(r"```(?:yaml|yml)?[^\n]*\n(.*?)```", text, flags=re.DOTALL):
        try:
            data = yaml.safe_load(block)
        except yaml.YAMLError as error:
            logger.warning("YAML block not readable: %s", error)
            continue
        if isinstance(data, dict) and isinstance(data.get("evidence"), dict):
            evidence = data["evidence"]
        elif isinstance(data, dict) and eln is None:
            eln = data
    return eln, evidence


def file_key(name):
    """A file name as llm_hub links it, without its `<hid>_` prefix and extension."""
    name = re.sub(r"^\d+_", "", Path(str(name).strip()).name)
    stem = name.rsplit(".", 1)[0] if "." in name else name
    return re.sub(r"[^\w\-]", "_", stem)


def valid_formula(formula):
    """True for element symbols with counts, ()/[] groups with counts and a hydrate dot (·)."""
    tokens = re.findall(r"[A-Z][a-z]?|\d+(?:\.\d+)?|[()\[\]·]|.", formula)
    stack, previous = [], None
    for token in tokens:
        if token[0].isupper():
            if token not in ELEMENTS:
                return False
        elif token[0].isdigit():
            # A count follows an element or a group; a leading count only starts a hydrate part.
            if previous not in (None, "·") and not (previous[0].isupper() or previous in ")]"):
                return False
        elif token in "([":
            stack.append(token)
        elif token in ")]":
            if not stack or {"(": ")", "[": "]"}[stack.pop()] != token or previous in (None, "(", "["):
                return False
        elif token == "·":
            if previous in (None, "·") or stack:
                return False
        else:
            return False
        previous = token
    return bool(tokens) and not stack and previous != "·" and any(t[0].isupper() for t in tokens)


def clean_formula(value, changes):
    formula = re.sub(r"\s+", "", value)
    if formula != value:
        changes.append(f"spaces removed: {value!r} -> {formula!r}")
    if not valid_formula(formula):
        raise Dropped("invalid_formula")
    return formula


def clean_name(value, changes):
    """ASCII name: Greek letters spelled out, dashes as `-`, a formula split by spaces joined."""
    name = value
    for dash in DASHES:
        name = name.replace(dash, "-")
    for letter, word in GREEK.items():
        name = re.sub(rf"{letter}-?", f"{word}-", name)
    # "Bi2 Pd" -> "Bi2Pd": the left word ends with an element and count, the right word is a formula.
    formula_word = r"(?:[A-Z][a-z]?\d*(?:\.\d+)?)+(?:\(\d+\))?"

    def join(match):
        left, right = match.group(1), match.group(2)
        symbols = re.findall(r"[A-Z][a-z]?", re.sub(r"\(\d+\)", "", right))
        return left + right if all(symbol in ELEMENTS for symbol in symbols) else match.group(0)
    name = re.sub(rf"([A-Z][a-z]?\d+)\s+({formula_word})(?=\s|$|,)", join, name)
    if name != value:
        changes.append(f"name made ASCII: {value!r} -> {name!r}")
    if not name.isascii():
        raise Dropped("non_ascii_name")
    return name


def allowed_value(value, allowed, changes):
    for option in allowed:
        if norm(value).lower() == option.lower():
            if value != option:
                changes.append(f"written as the allowed value: {value!r} -> {option!r}")
            return option
    raise Dropped("not_allowed")


def check_field(field, value, evidence, sources, rules, definition, examples):
    """The value to keep and its report entry; raises Dropped with the reason."""
    entry = {"value": value, "notes": [], "changes": []}
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        value = str(value)
    if not isinstance(value, str):
        raise Dropped("not_text")
    item = evidence.get(field)
    if not isinstance(item, dict) or not str(item.get("quote", "")).strip():
        raise Dropped("no_evidence")
    quote, named = norm(item["quote"]), file_key(item.get("file", ""))
    entry.update(file=item.get("file"), quote=item["quote"])
    if named in (file_key(INSTRUCTIONS_NAME), file_key(BASE_ELN_NAME)):
        raise Dropped("quote_from_instructions_or_base_eln")
    found = [source for source in sources if quote in norm(source["text"])]
    if not found:
        # PDF text can split a word ("Bi2 Pd"), so compare once more without any spaces.
        bare = quote.replace(" ", "")
        found = [source for source in sources if bare in norm(source["text"]).replace(" ", "")]
        if found:
            entry["notes"].append("quote found only when spaces are ignored")
    if not found:
        raise Dropped("quote_not_found")
    if named not in [file_key(source["name"]) for source in found]:
        entry["notes"].append(f"quote found in {found[0]['name']}, not in the named file")
    quote_text = norm(HEADER_KEY.sub(" ", quote))
    if not quote_text:
        raise Dropped("header_key_only")
    # Only the own quote, without header keys, counts: a short value like "3.3" is in most headers.
    if field in rules["literal"] and squash(value) not in squash(quote_text):
        raise Dropped("value_not_in_quote")
    if field == FORMULA_FIELD:
        value = clean_formula(value, entry["changes"])
    if field in NAME_FIELDS:
        value = clean_name(value, entry["changes"])
    if field in rules["exact"]:
        value = allowed_value(value, rules["definitions"][definition][field], entry["changes"])
    for example in examples:
        if squash(example) in squash(value) and not any(
            squash(example) in squash(source["text"]) for source in sources
        ):
            raise Dropped(f"copied_from_instructions: {example!r}")
    entry["value"] = value
    return value, entry


def merge(rule_eln, answer_text, bundle, raw_file, rules):
    """The merged ELN and the report."""
    definition = definition_of(rule_eln, rules)
    report = {
        "schema_version": 1,
        "sha256": bundle["raw"]["sha256"],
        "raw_path": bundle["raw"]["path"],
        "record_id": bundle["record"]["record_id"],
        "definition": definition,
        "status": "merged",
        "answer_sha256": None,
        "instructions": None,
        "kept": [],
        "dropped": [],
        "fields": {},
    }
    if answer_text is None:
        report["status"] = "no_answer"
        return copy.deepcopy(rule_eln), report
    report["answer_sha256"] = hashlib.sha256(answer_text.encode("utf-8")).hexdigest()
    if not definition:
        report["status"] = "no_definition"
        return copy.deepcopy(rule_eln), report
    instructions = load_instructions(definition, rules)
    report["instructions"] = {
        "path": rules["definitions"][definition]["instructions"],
        "sha256": hashlib.sha256(instructions.encode("utf-8")).hexdigest(),
    }
    answer_eln, evidence = parse_answer(answer_text)
    if answer_eln is None:
        report["status"] = "unparsable"
        return copy.deepcopy(rule_eln), report
    sources, _ = source_files(bundle, raw_file)
    examples = example_values(instructions)
    merged = copy.deepcopy(rule_eln)
    fields = [f for f in rules["fill"] if not has_path(rule_eln, f)] + rules["replace"]
    for field in fields:
        value = value_at(answer_eln, field)
        if value in (None, "", []):
            report["fields"][field] = {"status": "not_filled"}
            continue
        try:
            value, entry = check_field(field, value, evidence, sources, rules, definition, examples)
        except Dropped as reason:
            item = evidence.get(field) if isinstance(evidence.get(field), dict) else {}
            report["fields"][field] = {"status": "dropped", "reason": str(reason), "value": value,
                                       "file": item.get("file"), "quote": item.get("quote")}
            report["dropped"].append(field)
            continue
        set_path(merged, field, value)
        report["fields"][field] = {"status": "kept", **entry}
        report["kept"].append(field)
    return merged, report


def main():
    parser = argparse.ArgumentParser(description="Check an LLM ELN and merge it into the rule ELN.")
    parser.add_argument("bundle", help="Context bundle JSON of one raw file.")
    parser.add_argument("eln", help="ELN YAML of the rule builder.")
    parser.add_argument("--answer", help="The llm_hub answer (Markdown); without it the rule ELN is kept.")
    parser.add_argument("--answer-label", help="Name of the answer dataset (names the model), for the report.")
    parser.add_argument("--raw", help="The raw file; its .sxm header is a source, as for the LLM.")
    parser.add_argument("--rules", default=str(TOOL_DIR / "eln_rules.yaml"))
    parser.add_argument("--out-eln", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )
    bundle = json.loads(Path(args.bundle).read_text(encoding="utf-8"))
    rule_eln = yaml.safe_load(Path(args.eln).read_text(encoding="utf-8")) or {}
    rules = yaml.safe_load(Path(args.rules).read_text(encoding="utf-8"))
    answer = Path(args.answer).read_text(encoding="utf-8", errors="replace") if args.answer else None
    # A skipped llm_hub step (LLM off in the workflow) gives a dataset that holds only "null".
    if answer is not None and answer.strip() == "null":
        answer = None
    merged, report = merge(rule_eln, answer, bundle, args.raw, rules)
    report["answer_label"] = args.answer_label
    header = (
        f"# ELN for {report['raw_path']}: rule ELN + LLM values that passed the checks.\n"
        f"# LLM step: {report['status']}; kept: {', '.join(report['kept']) or 'none'}; "
        f"dropped: {', '.join(report['dropped']) or 'none'}.\n"
    )
    Path(args.out_eln).write_text(
        header + yaml.safe_dump(merged, allow_unicode=True, sort_keys=False), encoding="utf-8"
    )
    Path(args.report).write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("%s: %s; kept %s; dropped %s", report["raw_path"], report["status"],
                report["kept"], report["dropped"])


if __name__ == "__main__":
    main()
