import argparse
import pooch
import shutil
from pathlib import Path

def download_files(zenodo_doi_or_url, file_names, output_path, hash_dict=None):
    if zenodo_doi_or_url.startswith("https://"):
        zenodo_doi = zenodo_doi_or_url.split("doi.org/")[-1]
    else:
        zenodo_doi = zenodo_doi_or_url

    slug = zenodo_doi.replace("/", "_").replace(".", "_")
    root = Path(output_path) / slug
    root.mkdir(parents=True, exist_ok=True)

    if hash_dict:
        registry = {name: hash_dict[name] for name in file_names if name in hash_dict}

        pooch_obj = pooch.create(
            path=pooch.os_cache("pooch") / slug,
            base_url=f"doi:{zenodo_doi}",
            registry=registry,
            retry_if_failed=10,
            allow_updates=True,
        )
        fetch_files = lambda name: pooch_obj.fetch(name)
    else:
        fetch_files = lambda name: pooch.retrieve(
            url=f"https://zenodo.org/record/{zenodo_doi.split('.')[-1]}/files/{name}?download=1",
            known_hash=None,
            path=Path(output_path) / slug,
            fname=name
        )
    downloaded_files =[]
    for name in file_names:
        try:
            src = fetch_files(name)
            #shutil.copy(src, root / name)
            downloaded_files.append(str((root / name).resolve()))
        except Exception as e:
            print(e)


    return str(root.as_posix()),downloaded_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--zenodo_doi", required=True)
    parser.add_argument("--file_list", required=True)
    parser.add_argument("--hash_dict", required=False)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--output_path_file", required=True)

    args = parser.parse_args()

    with open(args.file_list) as f:
        file_names = [line.strip() for line in f if line.strip()]

    result_path,downloaded_files = download_files(args.zenodo_doi, file_names, args.output_dir,args.hash_dict)

    with open(args.output_path_file, "w") as f:
        f.write(result_path)

