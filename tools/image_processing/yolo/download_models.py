import argparse
import json
import os
import requests

def download_model(model_name, folder_name, json_file):
    with open(json_file, 'r') as f:
        model_dict = json.load(f)
    
    url = model_dict[model_name]
    filename = f"{model_name}.pt"
    
    destination_path = os.path.join(folder_name, filename)

    print(f"Downloading model '{model_name}' from {url} to {destination_path}...")

    response = requests.get(url)

    with open(destination_path, 'wb') as f:
        f.write(response.content)

    print("Download complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download a YOLO model .pt file.")
    parser.add_argument("-m","--model_name", help="Name of the YOLO model.")
    parser.add_argument("-f","--foldername", help="Folder to save the model.")
    parser.add_argument("-c", "--url_json", help="url file with model URLs.")

    args = parser.parse_args()

    print(args.model_name)
    download_model(args.model_name, args.foldername, args.url_json)
