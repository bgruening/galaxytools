import argparse
import csv
import torch
import numpy as np
from transformers import AutoTokenizer, AutoModel, AutoConfig
from visualize_attn import create_visualizations_for_selected_sequences

def parse_fasta(fasta_file):
    sequences = []
    header = None
    current_seq = []
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header is not None:
                    sequences.append((header, "".join(current_seq)))
                header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if header is not None:
            sequences.append((header, "".join(current_seq)))
    return sequences

def remove_special_tokens(attention_mask):
    mask = attention_mask.clone()
    mask[:, 0] = 0
    mask[torch.arange(mask.size(0)), mask.sum(dim=1)-1] = 0
    return mask

def mean_pooling(last_hidden_state, attention_mask):
    """Attention-mask aware mean pooling."""
    mask = attention_mask.unsqueeze(-1).type_as(last_hidden_state)
    summed = torch.sum(last_hidden_state * mask, dim=1)
    counts = torch.clamp(mask.sum(dim=1), min=1e-9)
    return summed / counts


def max_pooling(last_hidden_state, attention_mask):
    """Mask-aware max pooling."""
    mask = attention_mask.unsqueeze(-1)
    masked = last_hidden_state.masked_fill(mask == 0, float("-inf"))
    return torch.max(masked, dim=1).values


def embed_sequences(sequences, tokenizer, model, device, batch_size=16, 
                    pooling="mean", attn_vis=False, attn_ids=None):
    headers, seqs = zip(*sequences)
    embeddings = []
    all_attentions = {}
    all_tokens = {}
    
    print("Running batched inference...")
    for i in range(0, len(seqs), batch_size):
        batch = list(seqs[i:i + batch_size])
        batch_headers = headers[i:i + batch_size]
        inputs = tokenizer(batch, return_tensors="pt", padding=True, truncation=True)
        inputs = {k: v.to(device) for k, v in inputs.items()}

        with torch.inference_mode():
            outputs = model(**inputs, output_attentions=attn_vis)
            last_hidden = outputs.last_hidden_state

            if attn_vis and attn_ids is not None:
                attentions = outputs.attentions
                batch_attn = torch.stack(attentions).cpu()
                input_ids = inputs["input_ids"].cpu()
                mask = inputs["attention_mask"].cpu()
                for j, seq_id in enumerate(batch_headers):
                    if attn_ids is not None and seq_id not in attn_ids and "all" not in attn_ids:
                        continue
                    seq_len = mask[j].sum().item()
                    all_attentions[seq_id] = batch_attn[
                        :, j, :, :seq_len, :seq_len
                    ]
                    
                    token_ids = input_ids[j, :seq_len].tolist()
                    all_tokens[seq_id] = tokenizer.convert_ids_to_tokens(token_ids)
            
            pool_mask = remove_special_tokens(inputs["attention_mask"])

            if pooling == "cls":
                emb = last_hidden[:, 0, :]
            elif pooling == "mean":
                emb = mean_pooling(last_hidden, pool_mask)
            elif pooling == "max":
                emb = max_pooling(last_hidden, pool_mask)
            else:
                raise ValueError("Pooling must be 'cls', 'mean' or 'max'")
        embeddings.append(emb.cpu().numpy())
    embeddings = np.vstack(embeddings)

    if not attn_vis:
        all_attentions = None
        all_tokens = None
    return headers, embeddings, all_attentions, all_tokens


def save_to_tsv(headers, embeddings, output_file):
    dim = embeddings.shape[1]
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["id"] + [f"dim_{i}" for i in range(dim)])
        for h, emb in zip(headers, embeddings):
            writer.writerow([h] + emb.tolist())


def main():
    parser = argparse.ArgumentParser(description="DNABERT-2 FASTA embedding script")
    parser.add_argument("--model", type=str,
                        default="../DNABERT-2-117M",
                        help="Model name or path")
    parser.add_argument("--fasta", type=str, required=True,
                        help="Input FASTA file")
    parser.add_argument("--output", type=str, default="embeddings.tsv",
                        help="Output TSV file")
    parser.add_argument("--batch_size", type=int, default=16)
    parser.add_argument("--pooling", choices=["cls", "mean", "max"], default="cls")
    parser.add_argument("--attn_vis", action="store_true",
                        help="Enable attention visualization")
    parser.add_argument("--attn_ids", type=lambda x: x.split(","), default=None,
                        help="Comma-separated sequence IDs for attention saving")
    parser.add_argument("--viz_dir", type=str, default="attention_viz",
                        help="Directory to save visualization outputs")

    args = parser.parse_args()

    print("Loading model...")
    tokenizer = AutoTokenizer.from_pretrained(args.model, trust_remote_code=True)
    config = AutoConfig.from_pretrained(args.model, trust_remote_code=True)
    config.pad_token_id = tokenizer.pad_token_id
    model = AutoModel.from_pretrained(args.model, config=config, trust_remote_code=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    model.eval()

    print("Reading FASTA...")
    sequences = parse_fasta(args.fasta)
    print(f"Loaded {len(sequences)} sequences")

    print("Embedding...")
    headers, embeddings, attention, tokens = embed_sequences(
        sequences,
        tokenizer,
        model,
        device=device,
        batch_size=args.batch_size,
        pooling=args.pooling,
        attn_vis= args.attn_vis,
        attn_ids=args.attn_ids
    )

    print("Saving TSV...")
    save_to_tsv(headers, embeddings, args.output)
    
    if args.attn_vis and args.attn_ids is not None:
        print("Generating visualizations...")
        selected = headers if "all" in args.attn_ids else args.attn_ids
        create_visualizations_for_selected_sequences(selected, attention, tokens, args.viz_dir)
    
    print("Done.")
    print("Output shape:", embeddings.shape)


if __name__ == "__main__":
    main()