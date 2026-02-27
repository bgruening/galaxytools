import gc
import os

import torchaudio
from pyannote.audio import Pipeline
from transformers import Wav2Vec2ForCTC, Wav2Vec2Processor
from whisperx import asr

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
}

model_dir = "whisperx_models"
for model_name in DEFAULT_ALIGN_MODELS_TORCH.values():
    try:
        model = getattr(torchaudio.pipelines, model_name).get_model(
            dl_kwargs={"model_dir": model_dir + "/torch/hub/checkpoints"}
        )
        del model
        print(f"Loaded model: {model_name}")
    except Exception as e:
        print(e)
    gc.collect()

for model_name in DEFAULT_ALIGN_MODELS_HF.values():
    try:
        model_a = Wav2Vec2Processor.from_pretrained(
            model_name, cache_dir=model_dir + "/huggingface/hub"
        )
        model_b = Wav2Vec2ForCTC.from_pretrained(
            model_name, cache_dir=model_dir + "/huggingface/hub"
        )
        print(f"Loaded model: {model_name}")
        del model_a
        del model_b
    except Exception as e:
        print(e)
    gc.collect()


_MODELS = {
    "tiny": "Systran/faster-whisper-tiny",
    "base": "Systran/faster-whisper-base",
    "small": "Systran/faster-whisper-small",
    "medium": "Systran/faster-whisper-medium",
    "large": "Systran/faster-whisper-large-v3",
    "turbo": "mobiuslabsgmbh/faster-whisper-large-v3-turbo",
}

for model_name in _MODELS.keys():
    try:
        model = asr.load_model(
            model_name,
            "cpu",
            compute_type="int8",
            download_root=model_dir + "/faster-whisper",
        )
        print(f"Downloaded model: {model_name}")
    except Exception as e:
        print(f"Failed to download model {model_name}: {e}")
    gc.collect()


# py-annotate
Pipeline.from_pretrained(
    "pyannote/speaker-diarization-3.1",
    cache_dir=model_dir + "/torch/pyannote",
    use_auth_token=os.getenv("HF_AUTH_TOKEN"),
)
