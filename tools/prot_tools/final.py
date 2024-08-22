import argparse
import ast
import copy
import gc
import os
import random
import re
import json
import plotly
import plotly.graph_objs as go
import matplotlib.pyplot as plt
from tqdm import tqdm

import numpy as np
import optuna
# from smac import HyperparameterOptimizationFacade as HPOFacade
# from smac import Scenario
# from smac.intensifier.hyperband import Hyperband
# from smac.multi_objective.parego import ParEGO
import pandas as pd
import torch
import torch.nn.functional as F
from Bio import SeqIO
from datasets import Dataset
from optuna.samplers import TPESampler
from sklearn.metrics import confusion_matrix, matthews_corrcoef, roc_auc_score, accuracy_score
from sklearn.model_selection import train_test_split
from torch import nn
from torch.nn import BCEWithLogitsLoss, CrossEntropyLoss, MSELoss
from torch.utils.data import DataLoader
from transformers import (
    T5EncoderModel,
    T5Tokenizer,
    Trainer,
    TrainerCallback,
    TrainerControl,
    TrainerState,
    TrainingArguments,
    set_seed,
)
# from ConfigSpace import ConfigurationSpace, Configuration
# from ConfigSpace.hyperparameters import UniformFloatHyperparameter, CategoricalHyperparameter, UniformIntegerHyperparameter
from transformers.models.t5.modeling_t5 import T5Config, T5PreTrainedModel, T5Stack
from transformers.modeling_outputs import SequenceClassifierOutput
from transformers.utils.model_parallel_utils import assert_device_map, get_device_map
# import warnings
# warnings.filterwarnings('ignore')


class MultiObjectiveEarlyStoppingAndSaveCallback(TrainerCallback):
    def __init__(self, early_stopping_patience, min_delta=0.001, output_dir='./model_output', filename='finetuned_model'):
        self.early_stopping_patience = early_stopping_patience
        self.min_delta = min_delta
        self.best_val_loss = float('inf')
        self.best_val_accuracy = float('-inf')
        self.wait = 0
        self.output_dir = output_dir
        self.filename = filename
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def on_evaluate(self, args, state: TrainerState, control: TrainerControl, **kwargs):
        val_loss = kwargs['metrics']['eval_loss']
        val_accuracy = kwargs['metrics']['eval_accuracy']
        model = kwargs['model']

        loss_improved = (self.best_val_loss - val_loss) > self.min_delta
        accuracy_improved = (val_accuracy - self.best_val_accuracy) > self.min_delta

        if loss_improved or accuracy_improved:
            self.best_val_loss = min(self.best_val_loss, val_loss)
            self.best_val_accuracy = max(self.best_val_accuracy, val_accuracy)
            self.wait = 0
            # Save the model as the best so far
            self.save_finetuned_parameters(model, os.path.join(self.output_dir, self.filename))
            print(f"Saved improved model to {self.output_dir}/{self.filename}")
        else:
            self.wait += 1
            if self.wait >= self.early_stopping_patience:
                control.should_training_stop = True
                print(f"Stopping early at epoch {state.epoch}: No improvement in loss or accuracy for {self.early_stopping_patience} evaluations.")
                
    def save_finetuned_parameters(self, model, filepath):
        # Create a dictionary to hold the non-frozen parameters
        non_frozen_params = {n: p for n, p in model.named_parameters() if p.requires_grad}
        # Save only the finetuned parameters 
        torch.save(non_frozen_params, filepath)

def load_data(file_path, file_format):
    print(f"Loading data from {file_path} in {file_format} format")
    
    if file_format == 'fasta':
        sequences = []
        for record in SeqIO.parse(file_path, "fasta"):
            description_parts = record.description.split("%")
            label = int(description_parts[-1].split("LABEL=")[1])
            sequences.append([record.name, str(record.seq), label])
        df = pd.DataFrame(sequences, columns=["name", "sequence", "label"])
    elif file_format == 'csv':
        df = pd.read_csv(file_path)
    else:
        raise ValueError(f"Unsupported file format: {file_format}")
    
    return df

def split_data(data, train_percent, test_percent=None, val_percent=None, seed=42):
    if test_percent is None:
        # Two file mode
        train_data, val_data = train_test_split(data, train_size=train_percent, stratify=data['label'], random_state=seed)
        return train_data, val_data
    else:
        # Single file mode
        train_val, test = train_test_split(data, test_size=test_percent, stratify=data['label'], random_state=seed)
        train_size = train_percent / (train_percent + val_percent)
        train, val = train_test_split(train_val, train_size=train_size, stratify=train_val['label'], random_state=seed)
        return train, val, test

def load_and_prepare_data(file_path, test_size=0.2, random_state=42):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        description_parts = record.description.split("%")
        label = int(description_parts[-1].split("LABEL=")[1])
        sequences.append([record.name, str(record.seq), label])

    df = pd.DataFrame(sequences, columns=["name", "sequence", "label"])
    
    train_df, valid_df = train_test_split(df, test_size=test_size, random_state=random_state)
    
    train_df = train_df[["sequence", "label"]]
    valid_df = valid_df[["sequence", "label"]]
    
    return train_df, valid_df

def preprocess_sequences(df):
    # Replace uncommon AAs with "X"
    df["sequence"] = df["sequence"].str.replace('|'.join(["O","B","U","Z"]), "X", regex=True)
    # Add spaces between each amino acid for PT5 to correctly use them
    df['sequence'] = df.apply(lambda row: " ".join(row["sequence"]), axis=1)
    return df

def create_dataset(tokenizer, seqs, labels):
    tokenized = tokenizer(seqs, max_length=1024, padding=True, truncation=True)
    dataset = Dataset.from_dict(tokenized)
    dataset = dataset.add_column("labels", labels)
    return dataset

def set_seeds(s):
    torch.manual_seed(s)
    np.random.seed(s)
    random.seed(s)
    set_seed(s)

def parse_optuna_param(param_str):
    try:
        if param_str.startswith('[') and param_str.endswith(']'):
            return ast.literal_eval(param_str)
        elif param_str.startswith('(') and param_str.endswith(')'):
            return tuple(map(float, param_str.strip('()').split(',')))
        else:
            raise ValueError(f"Invalid parameter format: {param_str}")
    except:
        raise ValueError(f"Unable to parse parameter: {param_str}")

def parse_param(param_str):
    print(f"Input to parse_param: {param_str}")
    # Handle Galaxy-specific encoding
    param_str = param_str.replace('__ob__', '[').replace('__cb__', ']')
    try:
        # Try to parse as JSON first (for normal values)
        return json.loads(param_str)
    except json.JSONDecodeError:
        # If not JSON, check if it's a range or list
        param_str = param_str.strip()
        if param_str.startswith('(') and param_str.endswith(')'):
            # It's a range
            values = [float(x.strip()) for x in param_str[1:-1].split(',')]
            if len(values) == 2:
                return tuple(values)  # (min, max)
            elif len(values) == 3:
                return tuple(values)  # (min, max, step)
            else:
                raise ValueError(f"Invalid range format: {param_str}")
        elif param_str.startswith('[') and param_str.endswith(']'):
            # It's a list of options
            return ast.literal_eval(param_str)
        else:
            # It's a single value
            try:
                return float(param_str)
            except ValueError:
                return param_str

def process_hyperparameters(args):
    if args.hyperparameter_method == 'manual':
        return {
            'lr': float(args.lr),
            'batch': int(args.batch),
            'accum': int(args.accum),
            'dropout': float(args.dropout),
            'weight_decay': float(args.weight_decay),
            'warmup_pct': float(args.warmup_pct),
            'epochs': int(args.epochs),
        }
    else:
        return {
            'lr': parse_param(args.lr),
            'batch': parse_param(args.batch),
            'accum': parse_param(args.accum),
            'dropout': parse_param(args.dropout),
            'weight_decay': parse_param(args.weight_decay),
            'warmup_pct': parse_param(args.warmup_pct),
            'epochs_per_trial': int(args.epochs_per_trial)
        }

def process_adaptation_method(args):
    if args.adaptation_method == 'lora':
        if args.hyperparameter_method == 'manual':
            return {
                'method': 'lora',
                'lora_rank': int(args.lora_rank),
                'lora_init_scale': float(args.lora_init_scale),
                'lora_scaling_rank': int(args.lora_scaling_rank)
            }
        else:
            return {
                'method': 'lora',
                'lora_rank': parse_param(args.lora_rank),
                'lora_init_scale': parse_param(args.lora_init_scale),
                'lora_scaling_rank': parse_param(args.lora_scaling_rank)
            }
    else:
        if args.hyperparameter_method == 'manual':
            return {
                'method': 'dora',
                'dora_rank': int(args.dora_rank),
                'dora_init_scale': float(args.dora_init_scale)
            }
        else:
            return {
                'method': 'dora',
                'dora_rank': parse_param(args.dora_rank),
                'dora_init_scale': parse_param(args.dora_init_scale)
            }

class LoRAConfig:
    def __init__(self, lora_rank=8, lora_init_scale=0.01, lora_scaling_rank=2):
        self.lora_rank = lora_rank
        self.lora_init_scale = lora_init_scale
        self.lora_modules = ".*SelfAttention|.*EncDecAttention"
        self.lora_layers = "q|k|v|o"
        self.trainable_param_names = ".*layer_norm.*|.*lora_[ab].*"
        self.lora_scaling_rank = lora_scaling_rank
        # lora_modules and lora_layers are specified with regular expressions
        # see https://www.w3schools.com/python/python_regex.asp for reference
        
class LoRALinear(nn.Module):
    def __init__(self, linear_layer, rank, scaling_rank, init_scale):
        super().__init__()
        self.in_features = linear_layer.in_features
        self.out_features = linear_layer.out_features
        self.rank = rank
        self.scaling_rank = scaling_rank
        self.weight = linear_layer.weight
        self.bias = linear_layer.bias
        if self.rank > 0:
            self.lora_a = nn.Parameter(torch.randn(rank, linear_layer.in_features) * init_scale)
            if init_scale < 0:
                self.lora_b = nn.Parameter(torch.randn(linear_layer.out_features, rank) * init_scale)
            else:
                self.lora_b = nn.Parameter(torch.zeros(linear_layer.out_features, rank))
        if self.scaling_rank:
            self.multi_lora_a = nn.Parameter(
                torch.ones(self.scaling_rank, linear_layer.in_features)
                + torch.randn(self.scaling_rank, linear_layer.in_features) * init_scale
            )
            if init_scale < 0:
                self.multi_lora_b = nn.Parameter(
                    torch.ones(linear_layer.out_features, self.scaling_rank)
                    + torch.randn(linear_layer.out_features, self.scaling_rank) * init_scale
                )
            else:
                self.multi_lora_b = nn.Parameter(torch.ones(linear_layer.out_features, self.scaling_rank))

    def forward(self, input):
        if self.scaling_rank == 1 and self.rank == 0:
            # parsimonious implementation for ia3 and lora scaling
            if self.multi_lora_a.requires_grad:
                hidden = F.linear((input * self.multi_lora_a.flatten()), self.weight, self.bias)
            else:
                hidden = F.linear(input, self.weight, self.bias)
            if self.multi_lora_b.requires_grad:
                hidden = hidden * self.multi_lora_b.flatten()
            return hidden
        else:
            # general implementation for lora (adding and scaling)
            weight = self.weight
            if self.scaling_rank:
                weight = weight * torch.matmul(self.multi_lora_b, self.multi_lora_a) / self.scaling_rank
            if self.rank:
                weight = weight + torch.matmul(self.lora_b, self.lora_a) / self.rank
            return F.linear(input, weight, self.bias)

    def extra_repr(self):
        return "in_features={}, out_features={}, bias={}, rank={}, scaling_rank={}".format(
            self.in_features, self.out_features, self.bias is not None, self.rank, self.scaling_rank
        )


def modify_with_lora(transformer, config):
    for m_name, module in dict(transformer.named_modules()).items():
        if re.fullmatch(config.lora_modules, m_name):
            for c_name, layer in dict(module.named_children()).items():
                if re.fullmatch(config.lora_layers, c_name):
                    assert isinstance(
                        layer, nn.Linear
                    ), f"LoRA can only be applied to torch.nn.Linear, but {layer} is {type(layer)}."
                    setattr(
                        module,
                        c_name,
                        LoRALinear(layer, config.lora_rank, config.lora_scaling_rank, config.lora_init_scale),
                    )
    return transformer

class DoRAConfig:
    def __init__(self, dora_rank=8, dora_init_scale=0.01):
        self.dora_rank = dora_rank
        self.dora_init_scale = dora_init_scale
        self.dora_modules = ".*SelfAttention|.*EncDecAttention"
        self.dora_layers = "q|k|v|o"
        self.trainable_param_names = ".*layer_norm.*|.*dora_diag.*"

class DoRALinear(nn.Module):
    def __init__(self, linear_layer, rank, init_scale):
        super().__init__()
        self.in_features = linear_layer.in_features
        self.out_features = linear_layer.out_features
        self.rank = rank
        self.weight = linear_layer.weight
        self.bias = linear_layer.bias
        if self.rank > 0:
            self.dora_diag = nn.Parameter(torch.randn(min(self.in_features, self.out_features)) * init_scale)

    def forward(self, input):
        weight = self.weight
        if self.rank:
            diag_matrix = torch.diag(self.dora_diag)  # Create a diagonal matrix from the parameters
            diag_matrix_padded = torch.zeros_like(weight)
            diag_matrix_padded[:diag_matrix.shape[0], :diag_matrix.shape[1]] = diag_matrix
            weight = weight + diag_matrix_padded  # Ensure correct size
        return F.linear(input, weight, self.bias)

    def extra_repr(self):
        return "in_features={}, out_features={}, bias={}, rank={}".format(
            self.in_features, self.out_features, self.bias is not None, self.rank
        )


def modify_with_dora(transformer, config):
    for m_name, module in dict(transformer.named_modules()).items():
        if re.fullmatch(config.dora_modules, m_name):
            for c_name, layer in dict(module.named_children()).items():
                if re.fullmatch(config.dora_layers, c_name):
                    assert isinstance(
                        layer, nn.Linear
                    ), f"DoRA can only be applied to torch.nn.Linear, but {layer} is {type(layer)}."
                    setattr(
                        module,
                        c_name,
                        DoRALinear(layer, config.dora_rank, config.dora_init_scale),
                    )
    return transformer

class ClassConfig:
    def __init__(self, dropout=0.7, num_labels=2):
        self.dropout_rate = dropout
        self.num_labels = num_labels

class T5EncoderClassificationHead(nn.Module):
    """Head for sentence-level classification tasks."""

    def __init__(self, config, class_config):
        super().__init__()
        self.dense = nn.Linear(config.hidden_size, config.hidden_size)
        self.dropout = nn.Dropout(class_config.dropout_rate)
        self.out_proj = nn.Linear(config.hidden_size, class_config.num_labels)
        
        # Trainable emphasis factor
        self.emphasis_factor = nn.Parameter(torch.tensor(1.0))
        
    def forward(self, hidden_states):
        seq_length = hidden_states.size(1)
        middle_idx = seq_length // 2
        middle_embedding = hidden_states[:, middle_idx, :]

        # Apply trainable emphasis factor
        emphasized_middle_embedding = middle_embedding * self.emphasis_factor

        # Combine with the average embedding
        average_embedding = torch.mean(hidden_states, dim=1)
        combined_embedding = emphasized_middle_embedding + average_embedding

        x = self.dropout(combined_embedding)
        x = self.dense(x)
        x = torch.tanh(x)
        x = self.dropout(x)
        logits = self.out_proj(x)
        return logits

class T5EncoderForSimpleSequenceClassification(T5PreTrainedModel):

    def __init__(self, config: T5Config, class_config):
        super().__init__(config)
        self.num_labels = class_config.num_labels
        self.config = config

        self.shared = nn.Embedding(config.vocab_size, config.d_model)

        encoder_config = copy.deepcopy(config)
        encoder_config.use_cache = False
        encoder_config.is_encoder_decoder = False
        self.encoder = T5Stack(encoder_config, self.shared)

        self.dropout = nn.Dropout(class_config.dropout_rate) 
        self.classifier = T5EncoderClassificationHead(config, class_config)

        # Initialize weights and apply final processing
        self.post_init()

        # Model parallel
        self.model_parallel = False
        self.device_map = None

    def parallelize(self, device_map=None):
        self.device_map = (
            get_device_map(len(self.encoder.block), range(torch.cuda.device_count()))
            if device_map is None
            else device_map
        )
        assert_device_map(self.device_map, len(self.encoder.block))
        self.encoder.parallelize(self.device_map)
        self.classifier = self.classifier.to(self.encoder.first_device)
        self.model_parallel = True

    def deparallelize(self):
        self.encoder.deparallelize()
        self.encoder = self.encoder.to("cpu")
        self.model_parallel = False
        self.device_map = None
        torch.cuda.empty_cache()

    def get_input_embeddings(self):
        return self.shared

    def set_input_embeddings(self, new_embeddings):
        self.shared = new_embeddings
        self.encoder.set_input_embeddings(new_embeddings)

    def get_encoder(self):
        return self.encoder

    def _prune_heads(self, heads_to_prune):
        """
        Prunes heads of the model. heads_to_prune: dict of {layer_num: list of heads to prune in this layer} See base
        class PreTrainedModel
        """
        for layer, heads in heads_to_prune.items():
            self.encoder.layer[layer].attention.prune_heads(heads)

    def forward(
        self,
        input_ids=None,
        attention_mask=None,
        head_mask=None,
        inputs_embeds=None,
        labels=None,
        output_attentions=None,
        output_hidden_states=None,
        return_dict=None,
    ):
        return_dict = return_dict if return_dict is not None else self.config.use_return_dict

        outputs = self.encoder(
            input_ids=input_ids,
            attention_mask=attention_mask,
            inputs_embeds=inputs_embeds,
            head_mask=head_mask,
            output_attentions=output_attentions,
            output_hidden_states=True,
            return_dict=return_dict,
        )

        hidden_states = outputs[0]
        logits = self.classifier(hidden_states)

        loss = None
        if labels is not None:
            if self.config.problem_type is None:
                if self.num_labels == 1:
                    self.config.problem_type = "regression"
                elif self.num_labels > 1 and (labels.dtype == torch.long or labels.dtype == torch.int):
                    self.config.problem_type = "single_label_classification"
                else:
                    self.config.problem_type = "multi_label_classification"

            if self.config.problem_type == "regression":
                loss_fct = MSELoss()
                if self.num_labels == 1:
                    loss = loss_fct(logits.squeeze(), labels.squeeze())
                else:
                    loss = loss_fct(logits, labels)
            elif self.config.problem_type == "single_label_classification":
                loss_fct = CrossEntropyLoss()
                loss = loss_fct(logits.view(-1, self.num_labels), labels.view(-1))
            elif self.config.problem_type == "multi_label_classification":
                loss_fct = BCEWithLogitsLoss()
                loss = loss_fct(logits, labels)
        if not return_dict:
            output = (logits,) + outputs[1:]
            return ((loss,) + output) if loss is not None else output

        return SequenceClassifierOutput(
            loss=loss,
            logits=logits,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions,
        )

def PT5_classification_model(num_labels, dropout, adaptation_method='lora', lora_rank=None, lora_init_scale=None, lora_scaling_rank=None, dora_rank=None, dora_init_scale=None):
    # Load PT5 and tokenizer
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50", cache_dir="/scratch/users/hai/huggingface_cache")
    tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", cache_dir="/scratch/users/hai/huggingface_cache") 
    
    # Create new Classifier model with PT5 dimensions
    class_config=ClassConfig(num_labels=num_labels, dropout=dropout)
    class_model=T5EncoderForSimpleSequenceClassification(model.config,class_config)
    
    # Set encoder and embedding weights to checkpoint weights
    class_model.shared=model.shared
    class_model.encoder=model.encoder    
    
    # Delete the checkpoint model
    model=class_model
    del class_model
    
    # Print number of trainable parameters
    model_parameters = filter(lambda p: p.requires_grad, model.parameters())
    params = sum([np.prod(p.size()) for p in model_parameters])
    print("ProtT5_Classfier\nTrainable Parameter: "+ str(params))    
 
    if adaptation_method == 'lora':
        print("USING LoRA")
        # Add model modification lora
        config = LoRAConfig(lora_rank=lora_rank, lora_init_scale=lora_init_scale, lora_scaling_rank=lora_scaling_rank)
        
        # Add LoRA layers
        model = modify_with_lora(model, config)
    elif adaptation_method == 'dora':
        print("USING DoRA")
        # Add model modification dora
        config = DoRAConfig(dora_rank=dora_rank, dora_init_scale=dora_init_scale)
        
        # Add DoRA layers
        model = modify_with_dora(model, config)
    else:
        raise ValueError(f"Unsupported adaptation method: {adaptation_method}")
    
    # Freeze Embeddings and Encoder (except LoRA)
    for (param_name, param) in model.shared.named_parameters():
                param.requires_grad = False
    for (param_name, param) in model.encoder.named_parameters():
                param.requires_grad = False       

    for (param_name, param) in model.named_parameters():
            if re.fullmatch(config.trainable_param_names, param_name):
                param.requires_grad = True

    # Print trainable Parameter          
    model_parameters = filter(lambda p: p.requires_grad, model.parameters())
    params = sum([np.prod(p.size()) for p in model_parameters])
    print("ProtT5_LoRA_Classfier\nTrainable Parameter: "+ str(params) + "\n")
    
    return model, tokenizer

def train_per_protein(
        train_dataset,
        valid_dataset,
        weight_decay,
        warmup_pct,
        num_labels=2,
        batch=4,
        accum=2,
        val_batch=16,
        epochs=1,
        lr=3e-4,
        seed=42,
        deepspeed=False,
        gpu=1,
        dropout=0.5,
        adaptation_method='lora',
        lora_rank=4,
        lora_init_scale=0.01,
        lora_scaling_rank=1,
        dora_rank=4,
        dora_init_scale=0.01
):
    # Set gpu device
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu-1)
    
    # Set all random seeds
    set_seeds(seed)
    
    # load model
    model, tokenizer = PT5_classification_model(
        num_labels=num_labels,
        dropout=dropout,
        adaptation_method=adaptation_method,
        lora_rank=lora_rank,
        lora_init_scale=lora_init_scale,
        lora_scaling_rank=lora_scaling_rank,
        dora_rank=dora_rank,
        dora_init_scale=dora_init_scale
    )
    # Huggingface Trainer arguments
    total_steps = epochs * len(train_dataset) // batch
    warmup_steps = int(warmup_pct * total_steps)
     
    # Define TrainingArguments
    args = TrainingArguments(
        output_dir='./results',
        evaluation_strategy='epoch',
        logging_strategy='epoch',
        save_strategy='no',
        learning_rate=lr,
        per_device_train_batch_size=batch,
        gradient_accumulation_steps=accum,
        num_train_epochs=epochs,
        weight_decay=weight_decay,
        warmup_steps=warmup_steps,
        load_best_model_at_end=False,
        seed=seed,
        push_to_hub=False,
        logging_dir='./logs',
    )

    def compute_metrics(eval_pred):
        predictions, labels = eval_pred.predictions, eval_pred.label_ids
        if isinstance(predictions, tuple):
            predictions = predictions[0]
        if predictions.ndim > 1 and predictions.shape[1] > 1:
            predictions = np.argmax(predictions, axis=1)
        accuracy = accuracy_score(labels, predictions)
        return {"accuracy": accuracy}
    
    trainer = Trainer(
        model,
        args,
        train_dataset=train_dataset,
        eval_dataset=valid_dataset,
        tokenizer=tokenizer,
        compute_metrics=compute_metrics,
        callbacks=[MultiObjectiveEarlyStoppingAndSaveCallback(
            early_stopping_patience=3,
            min_delta=0.001,
            output_dir='./model_output',
            filename='finetuned_model_all.pth'
        )],
    )    
        
    # Train model
    trainer.train()

    # Load the best model
    best_model_path = os.path.join('./model_output', 'finetuned_model_all.pth')
    if os.path.exists(best_model_path):
        state_dict = torch.load(best_model_path)
        model.load_state_dict(state_dict, strict=False)
        print(f"Loaded best model from {best_model_path}")
        
    # Evaluate the best model
    eval_results = trainer.evaluate()
    print(eval_results)
    
    torch.cuda.empty_cache()
    gc.collect()

    return tokenizer, model, trainer.state.log_history

def get_train_per_protein_history(history):
    # Get loss, val_loss, and the computed metric from history
    loss = [x['loss'] for x in history if 'loss' in x]
    val_loss = [x['eval_loss'] for x in history if 'eval_loss' in x]

    # Get spearman (for regression) or accuracy value (for classification)
    if [x['eval_spearmanr'] for x in history if 'eval_spearmanr' in x] != []:
        metric = [x['eval_spearmanr'] for x in history if 'eval_spearmanr' in x]
        metric_name = 'validation_spearmanr'
    else:
        metric = [x['eval_accuracy'] for x in history if 'eval_accuracy' in x]
        metric_name = 'validation_accuracy'

    epochs = [x['epoch'] for x in history if 'loss' in x]

    # Truncate all lists to the shortest length
    min_length = min(len(loss), len(val_loss), len(metric), len(epochs))
    loss = loss[:min_length]
    val_loss = val_loss[:min_length]
    metric = metric[:min_length]
    epochs = epochs[:min_length]

    # Create a figure with two y-axes
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax2 = ax1.twinx()

    # Plot loss and val_loss on the first y-axis
    line1 = ax1.plot(epochs, loss, label='train_loss')
    line2 = ax1.plot(epochs, val_loss, label='validation_loss')
    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('Loss')

    # Plot the computed metric on the second y-axis
    line3 = ax2.plot(epochs, metric, color='red', label=metric_name)
    ax2.set_ylabel('Metric')
    ax2.set_ylim([0, 1])

    # Add grid lines
    ax1.grid(True)
    ax2.grid(True)

    # Combine the lines from both y-axes and create a single legend
    lines = line1 + line2 + line3
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels, loc='lower left')

    # Show the plot
    plt.title("Training History for fine-tuning")
    plt.savefig(f"Training_History_with_metric.pdf")

    # Close the plot to free up memory
    plt.close(fig)

    return fig  # Return the figure object in case you want to use it elsewhere

def metrics_calculation(model, tokenizer, my_test):
    # Set the device to use
    device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
    model.to(device)

    # create Dataset
    test_set=create_dataset(tokenizer,list(my_test['sequence']),list(my_test['label']))
    # make compatible with torch DataLoader
    test_set = test_set.with_format("torch", device=device)

    # Create a dataloader for the test dataset
    test_dataloader = DataLoader(test_set, batch_size=16, shuffle=False)

    # Put the model in evaluation mode
    model.eval()

    # Make predictions on the test dataset
    raw_logits = []
    labels = []
    with torch.no_grad():
        for batch in tqdm(test_dataloader):
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch['attention_mask'].to(device)
            # add batch results (logits) to predictions
            raw_logits += model(input_ids, attention_mask=attention_mask).logits.tolist()
            labels += batch["labels"].tolist()

    # Convert logits to predictions
    raw_logits = np.array(raw_logits)
    predictions = np.argmax(raw_logits, axis=1)

    # Calculate metrics
    conf_matrix = confusion_matrix(labels, predictions)
    tn, fp, fn, tp = conf_matrix.ravel()

    mcc = matthews_corrcoef(labels, predictions)
    specificity = tn / (tn + fp)
    sensitivity = tp / (tp + fn)
    accuracy = accuracy_score(labels, predictions)
    roc_auc = roc_auc_score(labels, raw_logits[:, 1])  # Assuming binary classification, adjust accordingly


    # Create a pandas DataFrame for metrics
    metrics_df = pd.DataFrame({
        "Metric": ["MCC", "Specificity", "Sensitivity", "Accuracy", "ROC-AUC"],
        "Value": [mcc, specificity, sensitivity, accuracy, roc_auc]
    })

    metrics_table = "metrics_table.csv"
    metrics_df.to_csv(metrics_table, index=False)

def main():
    parser = argparse.ArgumentParser(description='Process and validate input data and hyperparameters')
    parser.add_argument('--format', choices=['csv', 'fasta'], required=True)
    parser.add_argument('--input', help='Input file for single file mode')
    parser.add_argument('--train_val', help='Train/val file for two file mode')
    parser.add_argument('--test', help='Test file for two file mode')
    parser.add_argument('--train_percent', type=float, required=True)
    parser.add_argument('--test_percent', type=float)
    parser.add_argument('--val_percent', type=float)
    parser.add_argument('--percentage_summary', required=True, help='Output file for percentage summary')
    
    parser.add_argument('--model', choices=['pt5-xl50', 'esm', 'pt5-bfd'], required=True)
    parser.add_argument('--hyperparameter_method', choices=['manual', 'auto'], required=True)
    parser.add_argument('--optimization_algorithm', choices=['optuna', 'smac3'])
    parser.add_argument('--n_trials', type=int)
    parser.add_argument('--lr', type=str, required=True)
    parser.add_argument('--batch', type=str, required=True)
    parser.add_argument('--accum', type=str, required=True)
    parser.add_argument('--dropout', type=str, required=True)
    parser.add_argument('--weight_decay', type=str, required=True)
    parser.add_argument('--warmup_pct', type=str, required=True)
    parser.add_argument('--epochs', type=int)
    parser.add_argument('--epochs_per_trial', type=int)

    parser.add_argument('--adaptation_method', choices=['lora', 'dora'], required=True)
    parser.add_argument('--lora_rank', type=str)
    parser.add_argument('--lora_init_scale', type=str)
    parser.add_argument('--lora_scaling_rank', type=str)
    parser.add_argument('--dora_rank', type=str)
    parser.add_argument('--dora_init_scale', type=str)

    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')

    args = parser.parse_args()

    # Set the random seed for reproducibility
    np.random.seed(args.seed)
    try:

        # Load and split data
        if args.input:
            data = load_data(args.input, args.format)
            train_data, val_data, test_data = split_data(data, args.train_percent, args.test_percent, args.val_percent, seed=args.seed)
            print(f"Single file mode: Train size: {len(train_data)}, Validation size: {len(val_data)}, Test size: {len(test_data)}")
        else:
            train_val_data = load_data(args.train_val, args.format)
            test_data = load_data(args.test, args.format)
            train_data, val_data = split_data(train_val_data, args.train_percent, seed=args.seed)
            print(f"Two file mode: Train size: {len(train_data)}, Validation size: {len(val_data)}, Test size: {len(test_data)}")

        # Print class distribution for each split
        for name, dataset in [("Train", train_data), ("Validation", val_data), ("Test", test_data)]:
            class_distribution = dataset['label'].value_counts(normalize=True)
            print(f"\n{name} set class distribution:")
            for label, proportion in class_distribution.items():
                print(f"  Class {label}: {proportion:.2%}")

        print(f"\nRandom seed used: {args.seed}")

        # Load and prepare data, for now we will use only the train/val split. The test_size is the validation size
        # train_df, valid_df = load_and_prepare_data(args.train_val, test_size=1-args.train_percent)
        train_df, valid_df = train_data, val_data
        test_df = test_data
        
        # Preprocess sequences
        train_df = preprocess_sequences(train_df)
        valid_df = preprocess_sequences(valid_df)
        test_df = preprocess_sequences(test_df)

        # Print information about the datasets
        print("\nDataset shapes after preprocessing:")
        print(f"Train: {train_df.shape}")
        print(f"Validation: {valid_df.shape}")
        print(f"Test: {test_df.shape}")

        print("\nFirst few rows of each dataset:")
        print("\nTrain:")
        print(train_df.head())
        print("\nValidation:")
        print(valid_df.head())
        print("\nTest:")
        print(test_df.head())

        # Initialize the tokenizer
        tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", cache_dir="/scratch/users/hai/huggingface_cache")

        # Create Datasets
        train_set = create_dataset(tokenizer, list(train_df['sequence']), list(train_df['label']))
        valid_set = create_dataset(tokenizer, list(valid_df['sequence']), list(valid_df['label']))
        
        
        hyperparams = process_hyperparameters(args)
        adaptation_params = process_adaptation_method(args)
        
        with open(args.percentage_summary, 'w') as f:
            f.write("\n".join([f"{k}: {v}" for k, v in hyperparams.items()]))
            f.write("\n")
            f.write("\n".join([f"{k}: {v}" for k, v in adaptation_params.items()]))
            f.write("\n")
        
        print("Hyperparameters:")
        print(hyperparams["batch"])
            
        # raise SystemExit(0)
        if args.hyperparameter_method == 'manual':
            # Train the model with manual hyperparameters
            tokenizer, model, history = train_per_protein(
                train_set,
                valid_set,
                num_labels=2,
                batch=hyperparams['batch'],
                accum=hyperparams['accum'],
                epochs=hyperparams['epochs'],
                seed=args.seed,
                lr=hyperparams['lr'],
                dropout=hyperparams['dropout'],
                weight_decay=hyperparams['weight_decay'],
                warmup_pct=hyperparams['warmup_pct'],
                adaptation_method=adaptation_params['method'],
                **{k: v for k, v in adaptation_params.items() if k != 'method'}
            )
            
            print("Training completed. Final results:")
            print(history[-1])
            
            get_train_per_protein_history(history)
            
            metrics_calculation(model, tokenizer, test_df)


        else:  # Automated hyperparameter optimization
            def objective(trial):
                # Hyperparameters to be optimized
                lr = trial.suggest_float('lr', hyperparams['lr'][0], hyperparams['lr'][1], log=True)
                batch = trial.suggest_categorical('batch', hyperparams['batch'])
                accum = trial.suggest_categorical('accum', hyperparams['accum'])
                dropout = trial.suggest_float('dropout', hyperparams['dropout'][0], hyperparams['dropout'][1])
                weight_decay = trial.suggest_float('weight_decay', hyperparams['weight_decay'][0], hyperparams['weight_decay'][1], log=True)
                warmup_pct = trial.suggest_float("warmup_pct", hyperparams['warmup_pct'][0], hyperparams['warmup_pct'][1])

                if adaptation_params['method'] == 'lora':
                    lora_rank = trial.suggest_int('lora_rank', adaptation_params['lora_rank'][0], adaptation_params['lora_rank'][1], step=adaptation_params['lora_rank'][2] if len(adaptation_params['lora_rank']) > 2 else 1)
                    lora_init_scale = trial.suggest_float('lora_init_scale', adaptation_params['lora_init_scale'][0], adaptation_params['lora_init_scale'][1], log=True)
                    lora_scaling_rank = trial.suggest_int('lora_scaling_rank', adaptation_params['lora_scaling_rank'][0], adaptation_params['lora_scaling_rank'][1])
                    adaptation_trial_params = {'method': adaptation_params['method'],'lora_rank': lora_rank, 'lora_init_scale': lora_init_scale, 'lora_scaling_rank': lora_scaling_rank}
                else:
                    dora_rank = trial.suggest_int('dora_rank', adaptation_params['dora_rank'][0], adaptation_params['dora_rank'][1], step=adaptation_params['dora_rank'][2] if len(adaptation_params['dora_rank']) > 2 else 1)
                    dora_init_scale = trial.suggest_float('dora_init_scale', adaptation_params['dora_init_scale'][0], adaptation_params['dora_init_scale'][1], log=True)
                    adaptation_trial_params = {'method': adaptation_params['method'], 'dora_rank': dora_rank, 'dora_init_scale': dora_init_scale}

                
                # Training and evaluation
                tokenizer, model, history = train_per_protein(
                    train_dataset=train_set, 
                    valid_dataset=valid_set, 
                    num_labels=2, 
                    batch=batch, 
                    accum=accum, 
                    epochs=hyperparams['epochs_per_trial'],
                    lr=lr,
                    dropout=dropout,
                    weight_decay=weight_decay,
                    warmup_pct=warmup_pct,
                    seed=args.seed,
                    adaptation_method=adaptation_trial_params['method'],
                    **{k: v for k, v in adaptation_trial_params.items() if k != 'method'}
                )
                
                # Extract the last validation accuracy and loss from the history
                val_accuracy = history[-1]['eval_accuracy']
                val_loss = history[-1]['eval_loss']
                return val_loss, val_accuracy

            if args.optimization_algorithm == 'optuna':
                # Create and run the Optuna study
                study = optuna.create_study(
                    directions=['minimize', 'maximize'],
                    storage="sqlite:///optuna_results.sqlite3",
                    study_name="protein_classification_optimization",
                    sampler=TPESampler(seed=args.seed)
                )
                study.optimize(objective, n_trials=args.n_trials)

                # Print Pareto front
                print("Pareto front:")
                for trial in study.best_trials:
                    print(f"Trial {trial.number}")
                    print(f"  Loss: {trial.values[0]}, Accuracy: {trial.values[1]}")
                    print("  Params:")
                    for key, value in trial.params.items():
                        print(f"    {key}: {value}")

                best_trial = study.best_trials[0]  # You might want to choose based on your criteria
                best_params = best_trial.params

            # elif args.optimization_algorithm == 'smac3':
                # # Define the configuration space for SMAC3
                # cs = ConfigurationSpace()
                # cs.add_hyperparameter(UniformFloatHyperparameter("lr", hyperparams['lr'][0], hyperparams['lr'][1], log=True))
                # cs.add_hyperparameter(CategoricalHyperparameter("batch", choices=hyperparams['batch']))
                # cs.add_hyperparameter(CategoricalHyperparameter("accum", choices=hyperparams['accum']))
                # cs.add_hyperparameter(UniformFloatHyperparameter("dropout", hyperparams['dropout'][0], hyperparams['dropout'][1]))
                # cs.add_hyperparameter(UniformFloatHyperparameter("weight_decay", hyperparams['weight_decay'][0], hyperparams['weight_decay'][1], log=True))
                # cs.add_hyperparameter(UniformFloatHyperparameter("warmup_pct", hyperparams['warmup_pct'][0], hyperparams['warmup_pct'][1]))

                # if adaptation_params['method'] == 'lora':
                #     cs.add_hyperparameter(UniformIntegerHyperparameter("lora_rank", adaptation_params['lora_rank'][0], adaptation_params['lora_rank'][1]))
                #     cs.add_hyperparameter(UniformFloatHyperparameter("lora_init_scale", adaptation_params['lora_init_scale'][0], adaptation_params['lora_init_scale'][1], log=True))
                #     cs.add_hyperparameter(UniformIntegerHyperparameter("lora_scaling_rank", adaptation_params['lora_scaling_rank'][0], adaptation_params['lora_scaling_rank'][1]))
                # else:
                #     cs.add_hyperparameter(UniformIntegerHyperparameter("dora_rank", adaptation_params['dora_rank'][0], adaptation_params['dora_rank'][1]))
                #     cs.add_hyperparameter(UniformFloatHyperparameter("dora_init_scale", adaptation_params['dora_init_scale'][0], adaptation_params['dora_init_scale'][1], log=True))

                # # SMAC scenario
                # scenario = Scenario({
                #     "run_obj": "quality",
                #     "runcount-limit": args.n_trials,
                #     "cs": cs,
                #     "deterministic": "true",
                #     "output_dir": "smac3_output",
                # })

                # def smac_objective(config):
                #     loss, accuracy = objective(config)
                #     return loss  # SMAC minimizes, so we return the loss

                # # Run SMAC optimization
                # smac = SMAC4HPO(scenario=scenario, rng=np.random.RandomState(args.seed),
                #                 tae_runner=smac_objective)
                # best_found_config = smac.optimize()

                # best_params = best_found_config.get_dictionary()
                # print("Best configuration found:")
                # print(best_params)

            print("\nTraining with best hyperparameters:")
            tokenizer, model, history = train_per_protein(
                train_set,
                valid_set,
                num_labels=2,
                batch=best_params['batch'],
                accum=best_params['accum'],
                epochs=args.epochs,
                seed=args.seed,
                lr=best_params['lr'],
                dropout=best_params['dropout'],
                weight_decay=best_params['weight_decay'],
                warmup_pct=best_params['warmup_pct'],
                **{k: v for k, v in best_params.items() if k in ['lora_rank', 'lora_init_scale', 'lora_scaling_rank', 'dora_rank', 'dora_init_scale']}
            )

            print("Final training completed. Results:")
            print(history[-1])
    except ValueError as e:
        print(f"Error: {str(e)}")
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
if __name__ == "__main__":
    main()