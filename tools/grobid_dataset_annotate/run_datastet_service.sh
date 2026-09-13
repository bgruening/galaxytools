#!/usr/bin/env bash
set -euo pipefail

endpoint=""
input=""
output=""
add_paragraph_context="0"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --endpoint) endpoint="$2"; shift 2 ;;
        --input) input="$2"; shift 2 ;;
        --output) output="$2"; shift 2 ;;
        --add-paragraph-context) add_paragraph_context="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

if [[ -z "$endpoint" || -z "$input" || -z "$output" ]]; then
    echo "Missing required argument" >&2
    exit 2
fi

work_dir="$PWD"
datastet_config="$work_dir/datastet-config.yml"
grobid_config="$work_dir/grobid.yaml"
grobid_tmp="$work_dir/grobid-tmp"
mkdir -p "$grobid_tmp"
cp /opt/grobid/datastet/resources/config/config.yml "$datastet_config"
cp /opt/grobid/grobid-home/config/grobid.yaml "$grobid_config"
sed -i "s|tmpPath: \"tmp/\"|tmpPath: \"$grobid_tmp/\"|" "$datastet_config"
sed -i "s|temp: \"tmp\"|temp: \"$grobid_tmp\"|" "$grobid_config"
export GROBID_CONFIG_PATH="$grobid_config"

cd /opt/grobid/datastet
java --add-opens java.base/java.lang=ALL-UNNAMED -Dorg.grobid.config="$grobid_config" -Dorg.grobid.home=/opt/grobid/grobid-home -jar build/libs/datastet-0.8.1-onejar.jar server "$datastet_config" > "$work_dir/datastet-service.log" 2>&1 &
server_pid=$!
cd "$work_dir"
trap 'kill "$server_pid" 2>/dev/null || true; wait "$server_pid" 2>/dev/null || true' EXIT

for _ in $(seq 1 240); do
    if curl -fsS http://localhost:8060/service/isalive >/dev/null 2>&1; then
        break
    fi
    if ! kill -0 "$server_pid" 2>/dev/null; then
        echo "DataStet service exited before becoming ready" >&2
        tail -n 200 datastet-service.log >&2 || true
        exit 1
    fi
    sleep 2
done

if ! curl -fsS http://localhost:8060/service/isalive >/dev/null 2>&1; then
    echo "DataStet service did not become ready" >&2
    tail -n 200 datastet-service.log >&2 || true
    exit 1
fi

curl_args=(
    --fail
    --show-error
    --silent
    --form "input=@${input}"
    --output "$output"
)

if [[ "$endpoint" == "annotateDatasetPDF" ]]; then
    curl_args+=(--form "addParagraphContext=${add_paragraph_context}")
fi

curl "${curl_args[@]}" "http://localhost:8060/service/${endpoint}"

if [[ ! -s "$output" ]]; then
    echo "DataStet returned an empty response" >&2
    tail -n 200 datastet-service.log >&2 || true
    exit 1
fi
