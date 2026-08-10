import os
from typing import List

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from bertviz import head_view, model_view

def sanitize_filename(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in s)


def save_layer_avg_attention_html(
    seq_id: str,
    tokens: List[str],
    attn_probs: np.ndarray,   # [L,H,S,S]
    out_html: str
):
    L, H, S, _ = attn_probs.shape
    fig = go.Figure()

    positions = list(range(S))
    ticktext = tokens

    for l in range(L):
        mat = attn_probs[l].mean(axis=0)  # [S,S]
        fig.add_trace(
            go.Heatmap(
                z=mat,
                x=positions,
                y=positions,
                colorscale="Viridis",
                visible=(l == 0),
                colorbar=dict(title="Avg Attn") if l == 0 else None,
                customdata=np.array(tokens)[None, :].repeat(S, axis=0),
                hovertemplate=(
                    "Query=%{y}<br>"
                    "Key=%{x}<br>"
                    "Attention=%{z:.4f}<extra></extra>"
                ),
            )
        )

    buttons = []
    for l in range(L):
        visible = [False] * L
        visible[l] = True
        buttons.append(
            dict(
                label=f"Layer {l}",
                method="update",
                args=[
                    {"visible": visible},
                    {"title": f"{seq_id} - Layer {l} (average over heads)"},
                ],
            )
        )

    fig.update_layout(
        title=f"{seq_id} - Layer 0 (average over heads)",
        xaxis=dict(
            title="Key token",
            tickmode="array",
            tickvals=positions,
            ticktext=ticktext,
        ),
        yaxis=dict(
            title="Query token",
            tickmode="array",
            tickvals=positions,
            ticktext=ticktext,
            autorange="reversed",
        ),
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                showactive=True,
                x=1.02,
                y=1.0,
                xanchor="left",
                yanchor="top",
            )
        ],
        margin=dict(l=100, r=250, t=80, b=120),
        width=1200,
        height=900,
    )

    fig.write_html(out_html, include_plotlyjs="cdn")



def compute_attention_statistics(attn_probs, exclude_cls_sep=True):
    """
    attn_probs: [L, H, S, S]
    Returns:
      centrality: [S]
      entropy: [S]
      layer_entropy: [L, S]
      cls_outgoing: [S]
      key_mask: [S] bool
    """
    eps = 1e-12
    attn = np.asarray(attn_probs, dtype=np.float64)   # [L,H,S,S]
    L, H, S, _ = attn.shape

    # keys to keep
    key_mask = np.ones(S, dtype=bool)
    if exclude_cls_sep and S >= 1:
        key_mask[0] = False          # CLS
    if exclude_cls_sep and S >= 2:
        key_mask[-1] = False         # SEP

    # average heads
    attn_h = attn.mean(axis=1)       # [L,S,S]

    # mask excluded keys + renormalize each query distribution
    attn_h = attn_h * key_mask[None, None, :]
    row_sum = attn_h.sum(axis=-1, keepdims=True)
    attn_h = attn_h / np.clip(row_sum, eps, None)

    # incoming attention, normalized
    centrality = attn_h.mean(axis=0).sum(axis=0)      # [S]
    centrality = centrality / np.clip(centrality.sum(), eps, None)

    # entropy over keys
    ent_lq = -(attn_h * np.log(attn_h + eps)).sum(axis=-1)  # [L,S]

    # normalize by number of valid keys
    K = int(key_mask.sum())
    layer_entropy = ent_lq / np.log(max(K, 2))
    entropy = layer_entropy.mean(axis=0)                     # [S]

    # CLS outgoing distribution over keys after exclusion
    cls_outgoing = attn[:, :, 0, :].mean(axis=(0, 1))        # [S]
    cls_outgoing = cls_outgoing * key_mask
    cls_outgoing = cls_outgoing / np.clip(cls_outgoing.sum(), eps, None)

    return centrality, entropy, layer_entropy, cls_outgoing, key_mask


def save_attention_statistics_html(seq_id, tokens, attn_probs, out_html, top_k=50):
    centrality, entropy, layer_entropy, cls_outgoing, key_mask = compute_attention_statistics(
        attn_probs, exclude_cls_sep=True
    )

    centrality = np.asarray(centrality).reshape(-1)
    entropy = np.asarray(entropy).reshape(-1)
    cls_outgoing = np.asarray(cls_outgoing).reshape(-1)
    layer_entropy = np.asarray(layer_entropy)
    key_mask = np.asarray(key_mask).astype(bool).reshape(-1)

    L, S = layer_entropy.shape
    assert len(tokens) == S, f"tokens length ({len(tokens)}) != S ({S})"
    assert key_mask.shape[0] == S, f"key_mask length ({key_mask.shape[0]}) != S ({S})"

    positions = np.arange(S)
    labels = np.array([f"{i}:{tok}" for i, tok in enumerate(tokens)], dtype=object)

    valid_idx = positions[key_mask]
    if len(valid_idx) == 0:
        raise ValueError("No valid tokens left after excluding CLS/SEP.")

    def topk_indices(values, candidate_idx, k):
        if len(candidate_idx) == 0:
            return np.array([], dtype=int)
        k_eff = min(int(k), len(candidate_idx))
        local = np.argsort(values[candidate_idx])[-k_eff:][::-1]
        return candidate_idx[local]

    def bottomk_indices(values, candidate_idx, k):
        if len(candidate_idx) == 0:
            return np.array([], dtype=int)
        k_eff = min(int(k), len(candidate_idx))
        local = np.argsort(values[candidate_idx])[:k_eff]
        return candidate_idx[local]

    many = len(valid_idx) > top_k

    if many:
        idx_cent = topk_indices(centrality, valid_idx, top_k)
        idx_ent_hi = topk_indices(entropy, valid_idx, top_k)
        idx_ent_lo = bottomk_indices(entropy, valid_idx, top_k)
        idx_cls = topk_indices(cls_outgoing, valid_idx, top_k)

        n_rows = 5
        subplot_titles = [
            f"Token attention centrality (top {top_k})",
            f"Token attention entropy (top {top_k} highest)",
            f"Token attention entropy (top {top_k} lowest)",
            f"CLS outgoing attention (top {top_k})",
            "Entropy evolution across layers",
        ]
        specs = [[{"type": "bar"}], [{"type": "bar"}], [{"type": "bar"}], [{"type": "bar"}], [{"type": "heatmap"}]]
        row_heatmap = 5
    else:
        idx_cent = valid_idx
        idx_ent_hi = valid_idx
        idx_ent_lo = None
        idx_cls = valid_idx

        n_rows = 4
        subplot_titles = [
            "Token attention centrality",
            "Token attention entropy",
            "CLS outgoing attention",
            "Entropy evolution across layers",
        ]
        specs = [[{"type": "bar"}], [{"type": "bar"}], [{"type": "bar"}], [{"type": "heatmap"}]]
        row_heatmap = 4

    fig = make_subplots(
        rows=n_rows,
        cols=1,
        subplot_titles=subplot_titles,
        vertical_spacing=0.08 if n_rows == 5 else 0.11,
        specs=specs
    )

    def compact_x(idx):
        x = np.arange(len(idx))
        ticktext = labels[idx].tolist()
        return x, ticktext

    # 1) Centrality
    x_cent, txt_cent = compact_x(idx_cent)
    fig.add_trace(
        go.Bar(
            x=x_cent,
            y=centrality[idx_cent],
            customdata=np.array(txt_cent, dtype=object),
            text=[f"{v:.3f}" for v in centrality[idx_cent]],
            hovertemplate="<b>%{customdata}</b><br>Centrality=%{y:.4f}<extra></extra>",
        ),
        row=1, col=1
    )
    fig.update_xaxes(
        tickmode="array", tickvals=x_cent.tolist(), ticktext=txt_cent, tickangle=-60,
        row=1, col=1
    )

    # 2) Entropy (high)
    x_hi, txt_hi = compact_x(idx_ent_hi)
    fig.add_trace(
        go.Bar(
            x=x_hi,
            y=entropy[idx_ent_hi],
            customdata=np.array(txt_hi, dtype=object),
            text=[f"{v:.3f}" for v in entropy[idx_ent_hi]],
            hovertemplate="<b>%{customdata}</b><br>Entropy=%{y:.4f}<extra></extra>",
        ),
        row=2, col=1
    )
    fig.update_xaxes(
        tickmode="array", tickvals=x_hi.tolist(), ticktext=txt_hi, tickangle=-60,
        row=2, col=1
    )

    current_row = 3

    # 3) Entropy (low) only for many tokens
    if idx_ent_lo is not None:
        x_lo, txt_lo = compact_x(idx_ent_lo)
        fig.add_trace(
            go.Bar(
                x=x_lo,
                y=entropy[idx_ent_lo],
                customdata=np.array(txt_lo, dtype=object),
                text=[f"{v:.3f}" for v in entropy[idx_ent_lo]],
                hovertemplate="<b>%{customdata}</b><br>Entropy=%{y:.4f}<extra></extra>",
            ),
            row=3, col=1
        )
        fig.update_xaxes(
            tickmode="array", tickvals=x_lo.tolist(), ticktext=txt_lo, tickangle=-60,
            row=3, col=1
        )
        current_row = 4

    # 4) CLS outgoing
    x_cls, txt_cls = compact_x(idx_cls)
    fig.add_trace(
        go.Bar(
            x=x_cls,
            y=cls_outgoing[idx_cls],
            customdata=np.array(txt_cls, dtype=object),
            text=[f"{v:.3f}" for v in cls_outgoing[idx_cls]],
            hovertemplate="<b>%{customdata}</b><br>CLS→token attn=%{y:.4f}<extra></extra>",
        ),
        row=current_row, col=1
    )
    fig.update_xaxes(
        tickmode="array", tickvals=x_cls.tolist(), ticktext=txt_cls, tickangle=-60,
        row=current_row, col=1
    )

    # 5) Heatmap
    y_layers = [f"Layer {i}" for i in range(L)]
    hm_z = layer_entropy[:, valid_idx]                      # [L, V]
    hm_labels = labels[valid_idx]                           # [V]
    hm_x = np.arange(len(valid_idx))

    fig.add_trace(
        go.Heatmap(
            z=hm_z,
            x=hm_x,
            y=y_layers,
            colorscale="Viridis",
            colorbar=dict(
                title="Entropy",
                len=1.0,
                y=0.5,
                yanchor="middle",
                thickness=18
            ),
            customdata=np.tile(hm_labels, (L, 1)),
            hovertemplate="Layer %{y}<br>Token %{customdata}<br>Entropy=%{z:.3f}<extra></extra>",
        ),
        row=row_heatmap, col=1
    )

    if len(hm_x) > 80:
        step = max(1, len(hm_x) // 40)
        hm_tickvals = hm_x[::step]
    else:
        hm_tickvals = hm_x

    fig.update_xaxes(
        tickmode="array",
        tickvals=hm_tickvals.tolist(),
        ticktext=hm_labels[hm_tickvals].tolist(),
        tickangle=-60,
        row=row_heatmap, col=1
    )

    fig.update_layout(
        title=f"{seq_id} - Attention statistics",
        showlegend=False,
        width=1450,
        height=2000 if n_rows == 5 else 1650,
        margin=dict(l=120, r=120, t=100, b=150),
    )

    fig.write_html(out_html, include_plotlyjs="cdn")


def save_attention_matrix_tsv(
    tokens: List[str],
    attn_probs: np.ndarray,   # [L,H,S,S]
    out_path: str
):
    """
    Saves averaged attention matrix (over layers + heads)
    as a TSV with tokens as row/column labels.
    """
    attn_probs = np.asarray(attn_probs)
    # average over layers and heads
    attn_avg = attn_probs.mean(axis=(0, 1))   # [S,S]

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("token\t" + "\t".join(tokens) + "\n")
        for i, tok in enumerate(tokens):
            row = "\t".join(f"{x:.6f}" for x in attn_avg[i])
            f.write(f"{tok}\t{row}\n")

def save_attention_matrix_sheet(
    tokens: List[str],
    attn_probs: np.ndarray,   # [L,H,S,S]
    out_path: str
):
    """
    Save attention matrices into an Excel file with:
      - one sheet per layer
      - each sheet = average over heads (S x S)
      - rows = query tokens
      - columns = key tokens
    """

    attn_probs = np.asarray(attn_probs)
    L, H, S, _ = attn_probs.shape

    def sheet_name(l):
        return f"Layer_{l}"[:31]

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for l in range(L):
            mat = attn_probs[l].mean(axis=0)
            df = pd.DataFrame(
                mat,
                index=tokens,
                columns=tokens
            )
            df.index.name = "query \\ key"
            df.to_excel(writer, sheet_name=sheet_name(l))

def create_visualizations_for_selected_sequences(
    selected_indices, all_attention, all_tokens, viz_dir: str):

    os.makedirs(viz_dir, exist_ok=True)

    for seq_idx in selected_indices:
        seq_id = sanitize_filename(seq_idx)
        attn_probs = all_attention[seq_idx]
        tokens = all_tokens[seq_idx]

        # 1) Interactive per-layer average attention
        layer_avg_html = os.path.join(viz_dir, f"vis_{seq_id}_layer-avg-attention.html")
        save_layer_avg_attention_html(seq_idx, tokens, attn_probs, layer_avg_html)

        # 2+3) BertViz head and model view   
        if len(tokens) < 50:   
            attention = tuple(layer.unsqueeze(0) for layer in attn_probs)
            hv = head_view(attention, tokens, html_action="return")
            with open(os.path.join(viz_dir, f"vis_{seq_id}_bertviz-head-view.html"), "w", encoding="utf-8") as f:
                f.write(hv.data)

            mv = model_view(attention, tokens, html_action="return")
            with open(os.path.join(viz_dir, f"vis_{seq_id}_bertviz-model-view.html"), "w", encoding="utf-8") as f:
                f.write(mv.data)

        # 4) Bar plots
        stats_html = os.path.join(viz_dir, f"vis_{seq_id}_attention-statistics.html")
        save_attention_statistics_html(seq_idx, tokens, attn_probs, stats_html)
        
        # 5) Layer-Averaged Attention
        sheet_path = os.path.join(viz_dir, f"vis_{seq_id}_attn-avrg.xlsx")
        save_attention_matrix_sheet(tokens, attn_probs, sheet_path)