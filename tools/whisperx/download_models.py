"""Download the models of the Galaxy WhisperX tool into one directory.

The directory works for WhisperX 3.4.2 and 3.8.x and is the `path` of the
whisperx rows in the huggingface data table.
"""

import argparse
import os
import sys
import urllib.request
import zipfile
from pathlib import Path

import huggingface_hub

DEFAULT_ALIGN_MODELS_TORCH = {
    "en": "WAV2VEC2_ASR_BASE_960H",
    "fr": "VOXPOPULI_ASR_BASE_10K_FR",
    "de": "VOXPOPULI_ASR_BASE_10K_DE",
    "es": "VOXPOPULI_ASR_BASE_10K_ES",
    "it": "VOXPOPULI_ASR_BASE_10K_IT",
}

DEFAULT_ALIGN_MODELS_HF = {
    "ja": "jonatasgrosman/wav2vec2-large-xlsr-53-japanese",
    "zh": "jonatasgrosman/wav2vec2-large-xlsr-53-chinese-zh-cn",
    "nl": "jonatasgrosman/wav2vec2-large-xlsr-53-dutch",
    "uk": "Yehor/wav2vec2-xls-r-300m-uk-with-small-lm",
    "pt": "jonatasgrosman/wav2vec2-large-xlsr-53-portuguese",
    "ar": "jonatasgrosman/wav2vec2-large-xlsr-53-arabic",
    "cs": "comodoro/wav2vec2-xls-r-300m-cs-250",
    "ru": "jonatasgrosman/wav2vec2-large-xlsr-53-russian",
    "pl": "jonatasgrosman/wav2vec2-large-xlsr-53-polish",
    "hu": "jonatasgrosman/wav2vec2-large-xlsr-53-hungarian",
    "id": "cahya/wav2vec2-large-xlsr-indonesian",
    "fi": "jonatasgrosman/wav2vec2-large-xlsr-53-finnish",
    "fa": "jonatasgrosman/wav2vec2-large-xlsr-53-persian",
    "el": "jonatasgrosman/wav2vec2-large-xlsr-53-greek",
    "tr": "mpoyraz/wav2vec2-xls-r-300m-cv7-turkish",
    "da": "saattrupdan/wav2vec2-xls-r-300m-ftspeech",
    "he": "imvladikon/wav2vec2-xls-r-300m-hebrew",
    "vi": "nguyenvulebinh/wav2vec2-base-vi",
    "ko": "kresnik/wav2vec2-large-xlsr-korean",
    "ur": "kingabzpro/wav2vec2-large-xls-r-300m-Urdu",
    "te": "anuragshas/wav2vec2-large-xlsr-53-telugu",
    "hi": "theainerd/Wav2Vec2-large-xlsr-hindi",
    "ca": "softcatala/wav2vec2-large-xlsr-catala",
    "ml": "gvs/wav2vec2-large-xlsr-malayalam",
    "no": "NbAiLab/nb-wav2vec2-1b-bokmaal-v2",
    "nn": "NbAiLab/nb-wav2vec2-1b-nynorsk",
    "sk": "comodoro/wav2vec2-xls-r-300m-sk-cv8",
    "sl": "anton-l/wav2vec2-large-xlsr-53-slovenian",
    "hr": "classla/wav2vec2-xls-r-parlaspeech-hr",
    "ro": "gigant/romanian-wav2vec2",
    "eu": "stefan-it/wav2vec2-large-xlsr-53-basque",
    "gl": "ifrz/wav2vec2-large-xlsr-galician",
    "ka": "xsway/wav2vec2-large-xlsr-georgian",
    "lv": "jimregan/wav2vec2-large-xlsr-latvian-cv",
    "tl": "Khalsuu/filipino-wav2vec2-l-xls-r-300m-official",
    "sv": "KBLab/wav2vec2-large-voxrex-swedish",
}

_MODELS = {
    "tiny": "Systran/faster-whisper-tiny",
    "base": "Systran/faster-whisper-base",
    "small": "Systran/faster-whisper-small",
    "medium": "Systran/faster-whisper-medium",
    "large": "Systran/faster-whisper-large-v3",
    "turbo": "mobiuslabsgmbh/faster-whisper-large-v3-turbo",
}

TORCHAUDIO_MODEL_FILES = {
    "WAV2VEC2_ASR_BASE_960H": "wav2vec2_fairseq_base_ls960_asr_ls960.pth",
    "VOXPOPULI_ASR_BASE_10K_FR": "wav2vec2_voxpopuli_base_10k_asr_fr.pt",
    "VOXPOPULI_ASR_BASE_10K_DE": "wav2vec2_voxpopuli_base_10k_asr_de.pt",
    "VOXPOPULI_ASR_BASE_10K_ES": "wav2vec2_voxpopuli_base_10k_asr_es.pt",
    "VOXPOPULI_ASR_BASE_10K_IT": "wav2vec2_voxpopuli_base_10k_asr_it.pt",
}

FAILED = []

TORCHAUDIO_BASE_URL = "https://download.pytorch.org/torchaudio/models/"


def download_torchaudio_models(model_dir):
    # 3.4.2 reads $TORCH_HOME/hub/checkpoints/, 3.8.x reads --model_dir (hub/)
    hub_dir = Path(model_dir) / "hub"
    checkpoints_dir = hub_dir / "checkpoints"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)
    for pipeline_name, filename in TORCHAUDIO_MODEL_FILES.items():
        dest = checkpoints_dir / filename
        link = hub_dir / filename
        if not link.exists():
            link.unlink(missing_ok=True)
            link.symlink_to(Path("checkpoints") / filename)
        if dest.exists():
            print(f"Already exists, skipping: {filename}")
            continue
        url = TORCHAUDIO_BASE_URL + filename
        print(f"Downloading {pipeline_name} from {url} ...")
        try:
            urllib.request.urlretrieve(url, dest)
            print(f"Downloaded: {filename}")
        except Exception as e:
            dest.unlink(missing_ok=True)  # a partial file would be skipped on rerun
            print(f"Failed to download {pipeline_name}: {e}")
            FAILED.append(pipeline_name)


def download_hf_alignment_models(model_dir, token):
    cache_dir = str(Path(model_dir) / "hub")
    for lang, repo_id in DEFAULT_ALIGN_MODELS_HF.items():
        print(f"Downloading HF alignment model [{lang}]: {repo_id} ...")
        try:
            huggingface_hub.snapshot_download(repo_id, cache_dir=cache_dir, token=token)
            print(f"Downloaded: {repo_id}")
        except Exception as e:
            print(f"Failed to download {repo_id}: {e}")
            FAILED.append(repo_id)


def download_asr_models(model_dir, token):
    cache_dir = str(Path(model_dir) / "hub")
    for model_name, repo_id in _MODELS.items():
        print(f"Downloading ASR model [{model_name}]: {repo_id} ...")
        try:
            huggingface_hub.snapshot_download(
                repo_id,
                cache_dir=cache_dir,
                token=token,
                allow_patterns=[
                    "config.json",
                    "preprocessor_config.json",
                    "model.bin",
                    "tokenizer.json",
                    "vocabulary.*",
                ],
            )
            print(f"Downloaded: {model_name}")
        except Exception as e:
            print(f"Failed to download {model_name}: {e}")
            FAILED.append(model_name)


PYANNOTE_REPOS = [
    "pyannote/speaker-diarization-community-1",  # 3.8.x
    "pyannote/speaker-diarization-3.1",  # 3.4.2, with the next two
    "pyannote/segmentation-3.0",
    "pyannote/wespeaker-voxceleb-resnet34-LM",
]


def download_pyannote_model(model_dir, token):
    cache_dir = str(Path(model_dir) / "hub")
    for repo_id in PYANNOTE_REPOS:
        print(f"Downloading {repo_id} ...")
        try:
            huggingface_hub.snapshot_download(repo_id, cache_dir=cache_dir, token=token)
            print(f"Downloaded: {repo_id}")
        except Exception as e:
            print(f"Failed to download {repo_id}: {e}")
            FAILED.append(repo_id)


PUNKT_TAB_URL = "https://raw.githubusercontent.com/nltk/nltk_data/gh-pages/packages/tokenizers/punkt_tab.zip"


def download_punkt_tab(model_dir):
    # 3.8.x splits sentences with NLTK punkt_tab (NLTK_DATA=<model_dir>/nltk_data)
    tokenizers_dir = Path(model_dir) / "nltk_data" / "tokenizers"
    if (tokenizers_dir / "punkt_tab").is_dir():
        print("Already exists, skipping: punkt_tab")
        return
    tokenizers_dir.mkdir(parents=True, exist_ok=True)
    archive = tokenizers_dir / "punkt_tab.zip"
    print(f"Downloading punkt_tab from {PUNKT_TAB_URL} ...")
    try:
        urllib.request.urlretrieve(PUNKT_TAB_URL, archive)
        with zipfile.ZipFile(archive) as z:
            z.extractall(tokenizers_dir)
        print("Downloaded: punkt_tab")
    except Exception as e:
        print(f"Failed to download punkt_tab: {e}")
        FAILED.append("punkt_tab")
    finally:
        archive.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(
        description="Download WhisperX models (lightweight, no torch required)"
    )
    parser.add_argument(
        "--model-dir",
        default="whisperx_models",
        help="Directory to store downloaded models",
    )
    parser.add_argument(
        "--hf-token",
        default=os.getenv("HF_AUTH_TOKEN"),
        help="HuggingFace token for gated models",
    )
    args = parser.parse_args()

    model_dir = args.model_dir
    token = args.hf_token

    download_torchaudio_models(model_dir)
    download_hf_alignment_models(model_dir, token)
    download_asr_models(model_dir, token)
    download_pyannote_model(model_dir, token)
    download_punkt_tab(model_dir)

    if FAILED:
        # remove the loc rows for these models, or fix the cause and rerun
        print(f"Failed to download {len(FAILED)} model(s): {', '.join(FAILED)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
