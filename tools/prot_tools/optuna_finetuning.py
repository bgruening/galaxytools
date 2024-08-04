import os.path

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import BCEWithLogitsLoss, CrossEntropyLoss, MSELoss
from torch.utils.data import DataLoader

import re
import numpy as np
import pandas as pd
import copy

import transformers, datasets
from transformers.modeling_outputs import SequenceClassifierOutput
from transformers.models.t5.modeling_t5 import T5Config, T5PreTrainedModel, T5Stack
from transformers.utils.model_parallel_utils import assert_device_map, get_device_map
from transformers import T5EncoderModel, T5Tokenizer
from transformers import TrainingArguments, Trainer, set_seed

from evaluate import load
from datasets import Dataset

from tqdm import tqdm
import random

from scipy import stats
from sklearn.metrics import accuracy_score

import matplotlib.pyplot as plt
#!pip install umap-learn
import umap

from Bio import SeqIO
import pandas as pd

# sequences = []

# local_fasta_path = '../src/input_datasets/train_Pos_Neg_ST.fasta'

# # Load FASTA file using Biopython
# for record in SeqIO.parse(local_fasta_path, "fasta"):
#     # Split the description to extract label
#     description_parts = record.description.split("%")
#     label = int(description_parts[-1].split("LABEL=")[1])  # Extracting the numeric part of the label
#     sequences.append([record.name, str(record.seq), label])

# local_fasta_path = '../src/input_datasets/train_Pos_Neg_Y.fasta'

# for record in SeqIO.parse(local_fasta_path, "fasta"):
#     # Split the description to extract label
#     description_parts = record.description.split("%")
#     label = int(description_parts[-1].split("LABEL=")[1])  # Extracting the numeric part of the label
#     sequences.append([record.name, str(record.seq), label])

# # Create dataframe
# df = pd.DataFrame(sequences, columns=["name", "sequence", "label"])

# from sklearn.model_selection import train_test_split

# # Split the dataset into training and validation sets
# my_train, my_valid = train_test_split(df, test_size=0.2, random_state=42)

# my_train=my_train[["sequence", "label"]]
# my_valid=my_valid[["sequence","label"]]

# Modifies an existing transformer and introduce the LoRA layers

# class LoRAConfig:
#     def __init__(self, lora_rank=8, lora_init_scale=0.01, lora_scaling_rank=2):
#         self.lora_rank = lora_rank
#         self.lora_init_scale = lora_init_scale
#         self.lora_modules = ".*SelfAttention|.*EncDecAttention"
#         self.lora_layers = "q|k|v|o"
#         self.trainable_param_names = ".*layer_norm.*|.*lora_[ab].*"
#         self.lora_scaling_rank = lora_scaling_rank
#         # lora_modules and lora_layers are specified with regular expressions
#         # see https://www.w3schools.com/python/python_regex.asp for reference
        
# class LoRALinear(nn.Module):
#     def __init__(self, linear_layer, rank, scaling_rank, init_scale):
#         super().__init__()
#         self.in_features = linear_layer.in_features
#         self.out_features = linear_layer.out_features
#         self.rank = rank
#         self.scaling_rank = scaling_rank
#         self.weight = linear_layer.weight
#         self.bias = linear_layer.bias
#         if self.rank > 0:
#             self.lora_a = nn.Parameter(torch.randn(rank, linear_layer.in_features) * init_scale)
#             if init_scale < 0:
#                 self.lora_b = nn.Parameter(torch.randn(linear_layer.out_features, rank) * init_scale)
#             else:
#                 self.lora_b = nn.Parameter(torch.zeros(linear_layer.out_features, rank))
#         if self.scaling_rank:
#             self.multi_lora_a = nn.Parameter(
#                 torch.ones(self.scaling_rank, linear_layer.in_features)
#                 + torch.randn(self.scaling_rank, linear_layer.in_features) * init_scale
#             )
#             if init_scale < 0:
#                 self.multi_lora_b = nn.Parameter(
#                     torch.ones(linear_layer.out_features, self.scaling_rank)
#                     + torch.randn(linear_layer.out_features, self.scaling_rank) * init_scale
#                 )
#             else:
#                 self.multi_lora_b = nn.Parameter(torch.ones(linear_layer.out_features, self.scaling_rank))

#     def forward(self, input):
#         if self.scaling_rank == 1 and self.rank == 0:
#             # parsimonious implementation for ia3 and lora scaling
#             if self.multi_lora_a.requires_grad:
#                 hidden = F.linear((input * self.multi_lora_a.flatten()), self.weight, self.bias)
#             else:
#                 hidden = F.linear(input, self.weight, self.bias)
#             if self.multi_lora_b.requires_grad:
#                 hidden = hidden * self.multi_lora_b.flatten()
#             return hidden
#         else:
#             # general implementation for lora (adding and scaling)
#             weight = self.weight
#             if self.scaling_rank:
#                 weight = weight * torch.matmul(self.multi_lora_b, self.multi_lora_a) / self.scaling_rank
#             if self.rank:
#                 weight = weight + torch.matmul(self.lora_b, self.lora_a) / self.rank
#             return F.linear(input, weight, self.bias)

#     def extra_repr(self):
#         return "in_features={}, out_features={}, bias={}, rank={}, scaling_rank={}".format(
#             self.in_features, self.out_features, self.bias is not None, self.rank, self.scaling_rank
#         )


# def modify_with_lora(transformer, config):
#     for m_name, module in dict(transformer.named_modules()).items():
#         if re.fullmatch(config.lora_modules, m_name):
#             for c_name, layer in dict(module.named_children()).items():
#                 if re.fullmatch(config.lora_layers, c_name):
#                     assert isinstance(
#                         layer, nn.Linear
#                     ), f"LoRA can only be applied to torch.nn.Linear, but {layer} is {type(layer)}."
#                     setattr(
#                         module,
#                         c_name,
#                         LoRALinear(layer, config.lora_rank, config.lora_scaling_rank, config.lora_init_scale),
#                     )
#     return transformer

# class ClassConfig:
#     def __init__(self, dropout=0.7, num_labels=2):
#         self.dropout_rate = dropout
#         self.num_labels = num_labels

# class T5EncoderClassificationHead(nn.Module):
#     """Head for sentence-level classification tasks."""

#     def __init__(self, config, class_config):
#         super().__init__()
#         self.dense = nn.Linear(config.hidden_size, config.hidden_size)
#         self.dropout = nn.Dropout(class_config.dropout_rate)
#         self.out_proj = nn.Linear(config.hidden_size, class_config.num_labels)
        
#         # Trainable emphasis factor
#         self.emphasis_factor = nn.Parameter(torch.tensor(1.0))
        
#     def forward(self, hidden_states):
#         seq_length = hidden_states.size(1)
#         middle_idx = seq_length // 2
#         middle_embedding = hidden_states[:, middle_idx, :]

#         # Apply trainable emphasis factor
#         emphasized_middle_embedding = middle_embedding * self.emphasis_factor

#         # Combine with the average embedding
#         average_embedding = torch.mean(hidden_states, dim=1)
#         combined_embedding = emphasized_middle_embedding + average_embedding

#         x = self.dropout(combined_embedding)
#         x = self.dense(x)
#         x = torch.tanh(x)
#         x = self.dropout(x)
#         logits = self.out_proj(x)
#         return logits

    # def forward(self, hidden_states):

    #     hidden_states =  torch.mean(hidden_states,dim=1)  # avg embedding

    #     hidden_states = self.dropout(hidden_states)
    #     hidden_states = self.dense(hidden_states)
    #     hidden_states = torch.tanh(hidden_states)
    #     hidden_states = self.dropout(hidden_states)
    #     hidden_states = self.out_proj(hidden_states)
    #     return hidden_states
    
    # def forward(self, hidden_states):
    #     # Original sequence length and middle index
    #     seq_length = hidden_states.size(1)
    #     middle_idx = seq_length // 2

    #     # Extract the middle embedding vector
    #     middle_embedding = hidden_states[:, middle_idx, :]

    #     # Amplify the influence of the middle embedding
    #     amplified_middle_embedding = middle_embedding * 2

    #     # Combine with average to retain context
    #     average_embedding = torch.mean(hidden_states, dim=1)
    #     combined_embedding = 0.5 * amplified_middle_embedding + 0.5 * average_embedding

    #     # Classification layers
    #     x = self.dropout(combined_embedding)
    #     x = self.dense(x)
    #     x = torch.tanh(x)
    #     x = self.dropout(x)
    #     logits = self.out_proj(x)
    #     return logits


# class T5EncoderForSimpleSequenceClassification(T5PreTrainedModel):

#     def __init__(self, config: T5Config, class_config):
#         super().__init__(config)
#         self.num_labels = class_config.num_labels
#         self.config = config

#         self.shared = nn.Embedding(config.vocab_size, config.d_model)

#         encoder_config = copy.deepcopy(config)
#         encoder_config.use_cache = False
#         encoder_config.is_encoder_decoder = False
#         self.encoder = T5Stack(encoder_config, self.shared)

#         self.dropout = nn.Dropout(class_config.dropout_rate) 
#         self.classifier = T5EncoderClassificationHead(config, class_config)

#         # Initialize weights and apply final processing
#         self.post_init()

#         # Model parallel
#         self.model_parallel = False
#         self.device_map = None

#     def parallelize(self, device_map=None):
#         self.device_map = (
#             get_device_map(len(self.encoder.block), range(torch.cuda.device_count()))
#             if device_map is None
#             else device_map
#         )
#         assert_device_map(self.device_map, len(self.encoder.block))
#         self.encoder.parallelize(self.device_map)
#         self.classifier = self.classifier.to(self.encoder.first_device)
#         self.model_parallel = True

#     def deparallelize(self):
#         self.encoder.deparallelize()
#         self.encoder = self.encoder.to("cpu")
#         self.model_parallel = False
#         self.device_map = None
#         torch.cuda.empty_cache()

#     def get_input_embeddings(self):
#         return self.shared

#     def set_input_embeddings(self, new_embeddings):
#         self.shared = new_embeddings
#         self.encoder.set_input_embeddings(new_embeddings)

#     def get_encoder(self):
#         return self.encoder

#     def _prune_heads(self, heads_to_prune):
#         """
#         Prunes heads of the model. heads_to_prune: dict of {layer_num: list of heads to prune in this layer} See base
#         class PreTrainedModel
#         """
#         for layer, heads in heads_to_prune.items():
#             self.encoder.layer[layer].attention.prune_heads(heads)

#     def forward(
#         self,
#         input_ids=None,
#         attention_mask=None,
#         head_mask=None,
#         inputs_embeds=None,
#         labels=None,
#         output_attentions=None,
#         output_hidden_states=None,
#         return_dict=None,
#     ):
#         return_dict = return_dict if return_dict is not None else self.config.use_return_dict

#         outputs = self.encoder(
#             input_ids=input_ids,
#             attention_mask=attention_mask,
#             inputs_embeds=inputs_embeds,
#             head_mask=head_mask,
#             output_attentions=output_attentions,
#             output_hidden_states=True,
#             return_dict=return_dict,
#         )

#         hidden_states = outputs[0]
#         logits = self.classifier(hidden_states)

#         loss = None
#         if labels is not None:
#             if self.config.problem_type is None:
#                 if self.num_labels == 1:
#                     self.config.problem_type = "regression"
#                 elif self.num_labels > 1 and (labels.dtype == torch.long or labels.dtype == torch.int):
#                     self.config.problem_type = "single_label_classification"
#                 else:
#                     self.config.problem_type = "multi_label_classification"

#             if self.config.problem_type == "regression":
#                 loss_fct = MSELoss()
#                 if self.num_labels == 1:
#                     loss = loss_fct(logits.squeeze(), labels.squeeze())
#                 else:
#                     loss = loss_fct(logits, labels)
#             elif self.config.problem_type == "single_label_classification":
#                 loss_fct = CrossEntropyLoss()
#                 loss = loss_fct(logits.view(-1, self.num_labels), labels.view(-1))
#             elif self.config.problem_type == "multi_label_classification":
#                 loss_fct = BCEWithLogitsLoss()
#                 loss = loss_fct(logits, labels)
#         if not return_dict:
#             output = (logits,) + outputs[1:]
#             return ((loss,) + output) if loss is not None else output

#         return SequenceClassifierOutput(
#             loss=loss,
#             logits=logits,
#             hidden_states=outputs.hidden_states,
#             attentions=outputs.attentions,
#         )
        
# def PT5_classification_model(num_labels, dropout, lora_rank, lora_init_scale, lora_scaling_rank):
#     # Load PT5 and tokenizer
#     model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
#     tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50") 
    
#     # Create new Classifier model with PT5 dimensions
#     class_config=ClassConfig(num_labels=num_labels, dropout=dropout)
#     class_model=T5EncoderForSimpleSequenceClassification(model.config,class_config)
    
#     # Set encoder and embedding weights to checkpoint weights
#     class_model.shared=model.shared
#     class_model.encoder=model.encoder    
    
#     # Delete the checkpoint model
#     model=class_model
#     del class_model
    
#     # Print number of trainable parameters
#     model_parameters = filter(lambda p: p.requires_grad, model.parameters())
#     params = sum([np.prod(p.size()) for p in model_parameters])
#     print("ProtT5_Classfier\nTrainable Parameter: "+ str(params))    
 
#     # Add model modification lora
#     config = LoRAConfig(lora_rank=lora_rank, lora_init_scale=lora_init_scale, lora_scaling_rank=lora_scaling_rank)
    
#     # Add LoRA layers
#     model = modify_with_lora(model, config)
    
#     # Freeze Embeddings and Encoder (except LoRA)
#     for (param_name, param) in model.shared.named_parameters():
#                 param.requires_grad = False
#     for (param_name, param) in model.encoder.named_parameters():
#                 param.requires_grad = False       

#     for (param_name, param) in model.named_parameters():
#             if re.fullmatch(config.trainable_param_names, param_name):
#                 param.requires_grad = True

#     # Print trainable Parameter          
#     model_parameters = filter(lambda p: p.requires_grad, model.parameters())
#     params = sum([np.prod(p.size()) for p in model_parameters])
#     print("ProtT5_LoRA_Classfier\nTrainable Parameter: "+ str(params) + "\n")
    
#     return model, tokenizer

# from transformers import TrainerCallback, TrainerState, TrainerControl

# class EarlyStoppingCallback(TrainerCallback):
#     """Custom early stopping callback that can monitor loss or accuracy."""
    
#     def __init__(self, metric_name='eval_loss', early_stopping_patience=3, minimize=True):
#         """
#         Args:
#             metric_name (str): Metric to monitor, default 'eval_loss'.
#             early_stopping_patience (int): Number of checks with no improvement after which training will be stopped.
#             minimize (bool): Set to True if the metric should be minimized, False if it should be maximized.
#         """
#         self.metric_name = metric_name
#         self.early_stopping_patience = early_stopping_patience
#         self.early_stopping_counter = 0
#         self.minimize = minimize
#         self.best_metric = float('inf') if minimize else float('-inf')
    
#     def on_evaluate(self, args, state, control, **kwargs):
#         current_metric = kwargs['metrics'][self.metric_name]
        
#         if (self.minimize and current_metric < self.best_metric) or (not self.minimize and current_metric > self.best_metric):
#             self.best_metric = current_metric
#             self.early_stopping_counter = 0
#         else:
#             self.early_stopping_counter += 1
        
#         if self.early_stopping_counter >= self.early_stopping_patience:
#             control.should_training_stop = True
#             print(f'Stopping early! No improvement in {self.metric_name} for {self.early_stopping_patience} evaluation steps.')


# class MultiObjectiveEarlyStoppingCallback(TrainerCallback):
#     def __init__(self, early_stopping_patience, min_delta=0.001):
#         self.early_stopping_patience = early_stopping_patience
#         self.min_delta = min_delta
#         self.best_val_loss = float('inf')
#         self.best_val_accuracy = float('-inf')
#         self.wait = 0

#     def on_evaluate(self, args, state: TrainerState, control: TrainerControl, **kwargs):
#         # Extract current validation loss and accuracy
#         val_loss = kwargs['metrics']['eval_loss']
#         val_accuracy = kwargs['metrics']['eval_accuracy']

#         # Check if current loss and accuracy improved significantly
#         loss_improved = (self.best_val_loss - val_loss) > self.min_delta
#         accuracy_improved = (val_accuracy - self.best_val_accuracy) > self.min_delta

#         if loss_improved or accuracy_improved:
#             # Update best scores and reset wait time
#             self.best_val_loss = min(self.best_val_loss, val_loss)
#             self.best_val_accuracy = max(self.best_val_accuracy, val_accuracy)
#             self.wait = 0
#         else:
#             # If no improvement, increment the wait counter
#             self.wait += 1
#             if self.wait >= self.early_stopping_patience:
#                 # If wait exceeds the patience, stop training
#                 control.should_training_stop = True
#                 print(f"Stopping early at epoch {state.epoch}: No improvement in loss or accuracy for {self.early_stopping_patience} evaluations.")
                
# class MultiObjectiveEarlyStoppingAndSaveCallback(TrainerCallback):
#     def __init__(self, early_stopping_patience, min_delta=0.001, output_dir='./model_output', filename='finetuned_model'):
#         self.early_stopping_patience = early_stopping_patience
#         self.min_delta = min_delta
#         self.best_val_loss = float('inf')
#         self.best_val_accuracy = float('-inf')
#         self.wait = 0
#         self.output_dir = output_dir
#         self.filename = filename
#         if not os.path.exists(output_dir):
#             os.makedirs(output_dir)

#     def on_evaluate(self, args, state: TrainerState, control: TrainerControl, **kwargs):
#         val_loss = kwargs['metrics']['eval_loss']
#         val_accuracy = kwargs['metrics']['eval_accuracy']
#         model = kwargs['model']

#         loss_improved = (self.best_val_loss - val_loss) > self.min_delta
#         accuracy_improved = (val_accuracy - self.best_val_accuracy) > self.min_delta

#         if loss_improved or accuracy_improved:
#             self.best_val_loss = min(self.best_val_loss, val_loss)
#             self.best_val_accuracy = max(self.best_val_accuracy, val_accuracy)
#             self.wait = 0
#             # Save the model as the best so far
#             self.save_finetuned_parameters(model, os.path.join(self.output_dir, self.filename))
#             print(f"Saved improved model to {self.output_dir}/{self.filename}")
#         else:
#             self.wait += 1
#             if self.wait >= self.early_stopping_patience:
#                 control.should_training_stop = True
#                 print(f"Stopping early at epoch {state.epoch}: No improvement in loss or accuracy for {self.early_stopping_patience} evaluations.")
                
#     def save_finetuned_parameters(self, model, filepath):
#         # Create a dictionary to hold the non-frozen parameters
#         non_frozen_params = {n: p for n, p in model.named_parameters() if p.requires_grad}
#         # Save only the finetuned parameters 
#         torch.save(non_frozen_params, filepath)
        
# from sklearn.metrics import accuracy_score
# #!pip install seaborn
# import seaborn as sns
# import gc

# Set random seeds for reproducibility of your trainings run
# def set_seeds(s):
#     torch.manual_seed(s)
#     np.random.seed(s)
#     random.seed(s)
#     set_seed(s)

# def apply_umap(embeddings, n_components=2, min_dist=0.01):
#     umap_model = umap.UMAP(n_components=n_components)
#     umap_embeddings = umap_model.fit_transform(embeddings)
#     return umap_embeddings

# def plot_umap(embeddings, labels):
#     data = {"UMAP1": embeddings[:, 0], "UMAP2": embeddings[:, 1], "Label": labels}
#     df = pd.DataFrame(data)
    
#     plt.figure(figsize=(8, 6))
#     sns.scatterplot(x="UMAP1", y="UMAP2", hue="Label", data=df, palette={0: "blue", 1: "magenta"}, s=50, alpha=0.9)
#     plt.title("UMAP Visualization of Embeddings")
#     plt.savefig("../Plots/UMAP_Visualization_of_Embeddings_new.pdf")
#     plt.show()
    
# Main training fuction
# def train_per_protein(
#         train_dataset,         #training data
#         valid_dataset,         #validation data      
#         weight_decay,
#         warmup_pct,
#         num_labels= 2,    #1 for regression, >1 for classification
    
#         # effective training batch size is batch * accum
#         # we recommend an effective batch size of 8 
#         batch= 4,         #for training
#         accum= 2,         #gradient accumulation
    
#         val_batch = 16,   #batch size for evaluation
#         epochs=1,       #training epochs
#         lr= 3e-4,         #recommended learning rate
#         seed= 42,         #random seed
#         deepspeed=False,  #if gpu is large enough disable deepspeed for training speedup
#         gpu= 1,
#         dropout=0.5, #dropout rate
#          #L2 weight regularization
#         lora_rank=4,      #lora rank
#         lora_init_scale=0.01, #lora scaling rank
#         lora_scaling_rank=1,       #lora a
#         ):         #gpu selection (1 for first gpu)

#     # Set gpu device
#     os.environ["CUDA_VISIBLE_DEVICES"]=str(gpu-1)
    
#     # Set all random seeds
#     set_seeds(seed)
    
#     # load model
#     model, tokenizer = PT5_classification_model(num_labels=num_labels, dropout=dropout, lora_rank=lora_rank, lora_init_scale=lora_init_scale, lora_scaling_rank=lora_scaling_rank)

#     # Huggingface Trainer arguments
#     total_steps = epochs * len(train_dataset) // batch
#     warmup_steps = int(warmup_pct * total_steps)
     
#     # Define TrainingArguments
#     args = TrainingArguments(
#         output_dir='./results',              # where to save the model
#         evaluation_strategy='epoch',         # evaluation is done at the end of each epoch
#         logging_strategy='epoch',
#         save_strategy='no',
#         learning_rate=lr,                    # initial learning rate
#         per_device_train_batch_size=batch,   # batch size per device
#         gradient_accumulation_steps=accum,   # gradient accumulation steps
#         num_train_epochs=epochs,             # number of epochs to train
#         weight_decay=weight_decay,           # L2 weight regularization
#         warmup_steps=warmup_steps,           # 10% of total steps
#         load_best_model_at_end=False,         # load the best model at the end of training
#         seed=seed,                           # random seed
#         push_to_hub=False,                   # if you want to push model to the hub (Hugging Face Model Hub)
#         logging_dir='./logs',
#     )
#     # metric_for_best_model='eval_loss|accuracy'

#     # Metric definition for validation data
#     def compute_metrics(eval_pred):
#         predictions, labels = eval_pred.predictions, eval_pred.label_ids
#         # Check if predictions have the expected shape
#         if isinstance(predictions, tuple):
#             predictions = predictions[0]
#         if predictions.ndim > 1 and predictions.shape[1] > 1:
#             predictions = np.argmax(predictions, axis=1)
#         # Now, compute the metric (e.g., accuracy)
#         accuracy = accuracy_score(labels, predictions)
        
#         # Return the metric(s) as a dictionary
#         return {"accuracy": accuracy}
    
#     # For minimizing loss
#     early_stopping_loss = EarlyStoppingCallback(metric_name='eval_loss', early_stopping_patience=3, minimize=True)

#     # For maximizing accuracy
#     early_stopping_accuracy = EarlyStoppingCallback(metric_name='eval_accuracy', early_stopping_patience=3, minimize=False)
#     # Trainer          
#     trainer = Trainer(
#         model,
#         args,
#         train_dataset=train_dataset,
#         eval_dataset=valid_dataset,
#         tokenizer=tokenizer,
#         compute_metrics=compute_metrics,
#         callbacks=[MultiObjectiveEarlyStoppingAndSaveCallback(
#             early_stopping_patience=3,
#             min_delta=0.001,
#             output_dir='./model_output',
#             filename='finetuned_model_all.pth'
#         )],
#     )    

#     def get_embeddings(model, tokenizer, sequences, batch_size=32, device="cuda"):
#         embeddings = []
#         model = model.to(device)
#         model.eval()
    
#         # Iterate over the sequences in batches
#         for i in range(0, len(sequences), batch_size):
#             # Extract a batch of sequences
#             batch = sequences[i:i + batch_size]
    
#             # Tokenize the batch using the specified tokenizer and convert to PyTorch tensors
#             inputs = tokenizer(batch, return_tensors="pt", padding=True, truncation=True, max_length=512).to(device)
    
#             with torch.no_grad():
#                 # Forward pass through the model to obtain outputs
#                 outputs = model(**inputs)
    
#             # Extract hidden states from the second-to-last layer (penultimate layer)
#             hidden_states = outputs.hidden_states[-2].detach().cpu().numpy()
    
#             # Take the embeddings from the second-to-last layer
#             embeddings_from_layer = hidden_states[:, 0, :]
    
#             # Extend the list with the generated embeddings
#             embeddings.extend(embeddings_from_layer)
    
#             print(f"Batch {i // batch_size + 1}, Second-to-Last Layer Embeddings Shape: {embeddings_from_layer.shape}")
    
#         return np.array(embeddings)

        
#     # Train model
#     trainer.train()

#     # Get the best model
#     # model = trainer.model
#     # Ensure the best model is loaded
#     best_model_path = os.path.join('./model_output', 'finetuned_model_all.pth')
#     if os.path.exists(best_model_path):
#         state_dict = torch.load(best_model_path)
#         model.load_state_dict(state_dict, strict=False)
#         print(f"Loaded best model from {best_model_path}")
        
#     # Evaluate the best model
#     eval_results = trainer.evaluate()
#     print(eval_results)
    
#     # Print the current learning rate
#     # current_lr = trainer.optimizer.param_groups[0]['lr']
#     # print(f"Current learning rate: {current_lr}")
    
#     # valid_sequences = list(valid_dataset['sequence'])
#     # valid_embeddings = get_embeddings(model, tokenizer, valid_sequences)

#     # # Apply UMAP for dimensionality reduction
#     # umap_embeddings = apply_umap(valid_embeddings)

#     # # Plot UMAP embeddings
#     # labels = list(valid_dataset['label'])
#     # plot_umap(umap_embeddings, labels)
    
#     torch.cuda.empty_cache()
#     gc.collect()

#     return tokenizer, model, trainer.state.log_history


# # Dataset creation
# def create_dataset(tokenizer,seqs,labels):
#     tokenized = tokenizer(seqs, max_length=1024, padding=True, truncation=True)
#     dataset = Dataset.from_dict(tokenized)
#     dataset = dataset.add_column("labels", labels)

#     return dataset

# # Initialize the tokenizer
# tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50") 

# train_df = my_train
# valid_df = my_valid

# # Preprocess inputs
# # Replace uncommon AAs with "X"
# train_df["sequence"]=train_df["sequence"].str.replace('|'.join(["O","B","U","Z"]),"X",regex=True)
# valid_df["sequence"]=valid_df["sequence"].str.replace('|'.join(["O","B","U","Z"]),"X",regex=True)
# # Add spaces between each amino acid for PT5 to correctly use them
# train_df['sequence']=train_df.apply(lambda row : " ".join(row["sequence"]), axis = 1)
# valid_df['sequence']=valid_df.apply(lambda row : " ".join(row["sequence"]), axis = 1)

# # Create Datasets
# train_set=create_dataset(tokenizer,list(train_df['sequence']),list(train_df['label']))
# valid_set=create_dataset(tokenizer,list(valid_df['sequence']),list(valid_df['label']))

tokenizer, model, history = train_per_protein(train_set, valid_set, num_labels=2, batch=1, accum=2, epochs=20, seed=42, lr=0.00044666038459356726, dropout=0.6001375640608175, weight_decay=9.882084078511897e-05, warmup_pct=0.133784608876732, lora_rank=16, lora_init_scale=0.011516737020968842, lora_scaling_rank=3)