from __future__ import annotations

import torch
import torch.nn as nn

from pipeline.denovo.alphabet import MAX_CONTEXT, PAD_ID, VOCAB_SIZE


class HelixPeekLM(nn.Module):
    def __init__(
        self,
        vocab_size: int = VOCAB_SIZE,
        d_model: int = 96,
        nhead: int = 4,
        num_layers: int = 4,
        dim_feedforward: int = 384,
        dropout: float = 0.2,
        max_len: int = MAX_CONTEXT,
    ):
        super().__init__()
        if d_model % nhead:
            raise ValueError("d_model %d is not divisible by nhead %d" % (d_model, nhead))
        self.vocab_size = vocab_size
        self.d_model = d_model
        self.max_len = max_len

        self.token_embedding = nn.Embedding(vocab_size, d_model, padding_idx=PAD_ID)
        self.positional_embedding = nn.Parameter(torch.empty(1, max_len, d_model))
        nn.init.normal_(self.positional_embedding, mean=0.0, std=0.02)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            batch_first=True,
            norm_first=True,
            activation="gelu",
        )
        self.transformer_encoder = nn.TransformerEncoder(
            encoder_layer,
            num_layers=num_layers,
            norm=nn.LayerNorm(d_model),
            enable_nested_tensor=False,
        )
        self.lm_head = nn.Linear(d_model, vocab_size)

        self.register_buffer(
            "causal",
            torch.triu(torch.ones(max_len, max_len, dtype=torch.bool), diagonal=1),
            persistent=False,
        )

    def forward(self, token_ids: torch.Tensor) -> torch.Tensor:
        _, sequence_length = token_ids.shape
        if sequence_length > self.max_len:
            raise ValueError(
                "sequence of %d tokens exceeds the %d positions this model holds"
                % (sequence_length, self.max_len)
            )

        hidden_states = (
            self.token_embedding(token_ids)
            + self.positional_embedding[:, :sequence_length, :]
        )
        encoder_output = self.transformer_encoder(
            src=hidden_states,
            mask=self.causal[:sequence_length, :sequence_length],
            is_causal=True,
        )
        return self.lm_head(encoder_output)


def parameter_count(model: nn.Module) -> int:
    return sum(p.numel() for p in model.parameters())
