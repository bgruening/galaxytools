import argparse
import json
import re
import sys
import time
import zlib
from time import sleep
from urllib.parse import (
    parse_qs,
    urlencode,
    urlparse,
)
from xml.etree import ElementTree

import requests
from requests.adapters import (
    HTTPAdapter,
    Retry,
)


BATCH_SIZE = 50000  # Limit at UniProt is 100k
POLLING_INTERVAL = 5
API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        raise


def submit_id_mapping(from_db, to_db, ids):
    print(f"{from_db} {to_db}")
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] in ["NEW", "RUNNING"]:
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format in ["tsv", "gff"]:
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format in ["tsv", "gff"]:
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url, first):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if len(results) > 1 and file_format == "tsv" and not first:
        results = results[1:]
    if file_format == "xml":
        return merge_xml_results(results)
    return results


# print(results)
# {'results': [{'from': 'P05067', 'to': 'CHEMBL2487'}], 'failedIds': ['P12345']}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="retrieve uniprot mapping")
    subparsers = parser.add_subparsers(dest="tool")

    mapping = subparsers.add_parser("map")
    mapping.add_argument("f", help="from")
    mapping.add_argument("t", help="to")
    mapping.add_argument(
        "inp",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="input file (default: stdin)",
    )
    mapping.add_argument(
        "out",
        nargs="?",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="output file (default: stdout)",
    )
    mapping.add_argument("--format", default="tab", help="output format")
    mapping.add_argument("--fields", default="", help="Fields")

    retrieve = subparsers.add_parser("retrieve")
    retrieve.add_argument(
        "inp",
        metavar="in",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="input file (default: stdin)",
    )
    retrieve.add_argument(
        "out",
        nargs="?",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="output file (default: stdout)",
    )
    retrieve.add_argument("-f", "--format", help="specify output format", default="txt")
    mapping = subparsers.add_parser("menu")

    args = parser.parse_args()

    # code for auto generating the from - to conditional
    if args.tool == "menu":
        from lxml import etree

        request = session.get("https://rest.uniprot.org/configure/idmapping/fields")
        check_response(request)
        fields = request.json()

        tos = dict()
        from_cond = etree.Element("conditional", name="from_cond")
        from_select = etree.SubElement(
            from_cond, "param", name="from", type="select", label="Source database:"
        )

        rules = dict()
        for rule in fields["rules"]:
            rules[rule["ruleId"]] = rule["tos"]

        for group in fields["groups"]:
            group_name = group["groupName"]
            group_name = group_name.replace("databases", "DBs")
            for item in group["items"]:
                if item["to"]:
                    tos[item["name"]] = f"{group_name} - {item['displayName']}"

        for group in fields["groups"]:
            group_name = group["groupName"]
            group_name = group_name.replace("databases", "DBs")
            for item in group["items"]:
                if not item["from"]:
                    continue
                option = etree.SubElement(from_select, "option", value=item["name"])
                option.text = f"{group_name} - {item['displayName']}"
                when = etree.SubElement(from_cond, "when", value=item["name"])

                to_select = etree.SubElement(
                    when, "param", name="to", type="select", label="Target database:"
                )
                ruleId = item["ruleId"]
                for to in rules[ruleId]:
                    option = etree.SubElement(to_select, "option", value=to)
                    option.text = tos[to]
        etree.indent(from_cond, space="    ")
        print(etree.tostring(from_cond, pretty_print=True, encoding="unicode"))
        sys.exit(0)

    # get the IDs from the file as sorted list
    # (sorted is convenient for testing)
    query = set()
    for line in args.inp:
        query.add(line.strip())
    query = list(query)
    results = []
    first = True  # if False the header is removed
    while len(query) > 0:
        batch = query[:BATCH_SIZE]
        query = query[BATCH_SIZE:]
        print(f"processing {len(batch)} left {len(query)}")
        if args.tool == "map":
            job_id = submit_id_mapping(from_db=args.f, to_db=args.t, ids=batch)
        elif args.tool == "retrieve":
            job_id = submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=batch)

        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            link = f"{link}?format={args.format}"
            if args.tool == "map" and args.fields:
                link += f"&fields={args.fields}"
            print(link)
            results.extend(get_id_mapping_results_search(link, first))
            first = False
        print(f"got {len(results)} results so far")
        if len(query):
            sleep(5)

    if not isinstance(results, str):
        if args.format in ["fasta", "txt"]:
            results = "".join(results)
        else:
            results = "\n".join(results)
    args.out.write(f"{results}\n")
