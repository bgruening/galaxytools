#!/usr/bin/env bash
set -euo pipefail

endpoint=""
upload_mode=""
input=""
output=""
consolidate_header="0"
consolidate_citations="0"
include_raw_citations="0"
include_raw_affiliations="0"
start_page="-1"
end_page="-1"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --endpoint) endpoint="$2"; shift 2 ;;
        --upload-mode) upload_mode="$2"; shift 2 ;;
        --input) input="$2"; shift 2 ;;
        --output) output="$2"; shift 2 ;;
        --consolidate-header) consolidate_header="$2"; shift 2 ;;
        --consolidate-citations) consolidate_citations="$2"; shift 2 ;;
        --include-raw-citations) include_raw_citations="$2"; shift 2 ;;
        --include-raw-affiliations) include_raw_affiliations="$2"; shift 2 ;;
        --start-page) start_page="$2"; shift 2 ;;
        --end-page) end_page="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

if [[ -z "$endpoint" || -z "$upload_mode" || -z "$input" || -z "$output" ]]; then
    echo "Missing required argument" >&2
    exit 2
fi

work_dir="$PWD"
grobid_config="$work_dir/grobid.yaml"
grobid_tmp="$work_dir/grobid-tmp"
mkdir -p "$grobid_tmp"
cp /opt/grobid/grobid-home/config/grobid.yaml "$grobid_config"
# Galaxy runs Docker jobs as the invoking user. The full GROBID image stores
# DeLFT LMDB resources under root-owned paths, so prefer the Wapiti models.
sed -i 's/engine: "delft"/engine: "wapiti"/g' "$grobid_config"
sed -i "s|temp: \"tmp\"|temp: \"$grobid_tmp\"|" "$grobid_config"
export GROBID_CONFIG_PATH="$grobid_config"
export JAVA_OPTS="${JAVA_OPTS:-} -Dorg.grobid.config=$grobid_config -Dorg.grobid.home=/opt/grobid/grobid-home"

cd /opt/grobid
./grobid-service/bin/grobid-service > "$work_dir/grobid-service.log" 2>&1 &
server_pid=$!
cd "$work_dir"
trap 'kill "$server_pid" 2>/dev/null || true; wait "$server_pid" 2>/dev/null || true' EXIT

for _ in $(seq 1 180); do
    if curl -fsS http://localhost:8070/api/isalive >/dev/null 2>&1; then
        break
    fi
    if ! kill -0 "$server_pid" 2>/dev/null; then
        echo "GROBID service exited before becoming ready" >&2
        tail -n 200 grobid-service.log >&2 || true
        exit 1
    fi
    sleep 2
done

if ! curl -fsS http://localhost:8070/api/isalive >/dev/null 2>&1; then
    echo "GROBID service did not become ready" >&2
    tail -n 200 grobid-service.log >&2 || true
    exit 1
fi

curl_args=(
    --fail
    --show-error
    --silent
    --header "Accept: application/xml"
    --output "$output"
)

if [[ "$upload_mode" == "multipart" ]]; then
    curl_args+=(--form "input=@${input}")
    curl_args+=(
        --form "consolidateHeader=${consolidate_header}"
        --form "consolidateCitations=${consolidate_citations}"
        --form "includeRawCitations=${include_raw_citations}"
        --form "includeRawAffiliations=${include_raw_affiliations}"
    )
else
    curl_args+=(--data-urlencode "input@${input}")
    curl_args+=(
        --data "consolidateCitations=${consolidate_citations}"
        --data "includeRawCitations=${include_raw_citations}"
    )
fi

if [[ "$upload_mode" == "multipart" && "$start_page" != "-1" ]]; then
    curl_args+=(--form "start=${start_page}")
fi
if [[ "$upload_mode" == "multipart" && "$end_page" != "-1" ]]; then
    curl_args+=(--form "end=${end_page}")
fi

curl "${curl_args[@]}" "http://localhost:8070/api/${endpoint}"

if [[ ! -s "$output" ]]; then
    echo "GROBID returned an empty response" >&2
    tail -n 200 grobid-service.log >&2 || true
    exit 1
fi
