"""The 47,361-parameter 704D bounded native-geometry residual architecture.

The same residual architecture is reused by the GeoFuse-CPI V2.0 ATOM3D
native expert, where its target is an out-of-fold residual to the frozen V2.0
2D task head. Historical V16/V17 reference implementations remain read-only
and are not redistributed.
"""
from torch import nn
import torch


class ExplicitPairResidual(nn.Module):
    """Scaled, bounded correction to an independently calibrated 2D affinity prediction."""

    def __init__(self, n_dim=704):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_dim, 64),
            nn.LayerNorm(64),
            nn.SiLU(),
            nn.Linear(64, 32),
            nn.SiLU(),
            nn.Linear(32, 1),
        )

    def forward(self, x):
        return 0.5 * torch.tanh(self.net(x).squeeze(-1))


def make_calibrator(seed):
    """V17 LightGBM calibration settings (original old 349-source fourfold experiment)."""
    from lightgbm import LGBMRegressor

    return LGBMRegressor(
        n_estimators=175,
        learning_rate=0.04,
        num_leaves=9,
        max_depth=5,
        min_child_samples=12,
        reg_alpha=0.5,
        reg_lambda=8.0,
        colsample_bytree=0.85,
        subsample=0.9,
        subsample_freq=1,
        n_jobs=3,
        verbosity=-1,
        random_state=int(seed),
    )


def check_architecture():
    """Check byte-level architectural parity in dimensionality and parameter count."""
    net = ExplicitPairResidual()
    count = sum(p.numel() for p in net.parameters())
    if count != 47361:
        raise AssertionError(f"Unexpected residual model parameters: {count}")
    with torch.no_grad():
        out = net(torch.zeros((3, 704), dtype=torch.float32))
    if out.shape != (3,) or not torch.isfinite(out).all():
        raise AssertionError("Residual architecture forward check failed")
    return count
