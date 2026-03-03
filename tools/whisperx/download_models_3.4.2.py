import argparse
import os
import urllib.request
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

TORCHAUDIO_BASE_URL = "https://download.pytorch.org/torchaudio/models/"


def download_torchaudio_models(model_dir):
    checkpoints_dir = Path(model_dir) / "hub" / "checkpoints"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)
    for pipeline_name, filename in TORCHAUDIO_MODEL_FILES.items():
        dest = checkpoints_dir / filename
        if dest.exists():
            print(f"Already exists, skipping: {filename}")
            continue
        url = TORCHAUDIO_BASE_URL + filename
        print(f"Downloading {pipeline_name} from {url} ...")
        try:
            urllib.request.urlretrieve(url, dest)
            print(f"Downloaded: {filename}")
        except Exception as e:
            print(f"Failed to download {pipeline_name}: {e}")


def download_hf_alignment_models(model_dir, token):
    cache_dir = str(Path(model_dir) / "hub")
    for lang, repo_id in DEFAULT_ALIGN_MODELS_HF.items():
        print(f"Downloading HF alignment model [{lang}]: {repo_id} ...")
        try:
            huggingface_hub.snapshot_download(repo_id, cache_dir=cache_dir, token=token)
            print(f"Downloaded: {repo_id}")
        except Exception as e:
            print(f"Failed to download {repo_id}: {e}")


def download_asr_models(model_dir, token):
    cache_dir = str(Path(model_dir) / "hub")
    for model_name, repo_id in _MODELS.items():
        print(f"Downloading ASR model [{model_name}]: {repo_id} ...")
        try:
            huggingface_hub.snapshot_download(
                repo_id,
                cache_dir=cache_dir,
                token=token,
                allow_patterns=["config.json", "preprocessor_config.json", "model.bin", "tokenizer.json", "vocabulary.*"],
            )
            print(f"Downloaded: {model_name}")
        except Exception as e:
            print(f"Failed to download {model_name}: {e}")


PYANNOTE_REPOS = [
    "pyannote/speaker-diarization-3.1",
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


def main():
    parser = argparse.ArgumentParser(description="Download WhisperX models (lightweight, no torch required)")
    parser.add_argument("--model-dir", default="whisperx_models", help="Directory to store downloaded models")
    parser.add_argument("--hf-token", default=os.getenv("HF_AUTH_TOKEN"), help="HuggingFace token for gated models")
    args = parser.parse_args()

    model_dir = args.model_dir
    token = args.hf_token

    download_torchaudio_models(model_dir)
    download_hf_alignment_models(model_dir, token)
    download_asr_models(model_dir, token)
    download_pyannote_model(model_dir, token)


if __name__ == "__main__":
    main()
