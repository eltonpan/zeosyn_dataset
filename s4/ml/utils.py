"""Utility functions for ML purposes."""
from typing import Optional

import numpy
from matplotlib.colors import to_rgb


def weights2colors(weights: Optional[numpy.ndarray] = None, base_color: str = 'tab:blue'):
    """Convert weights to transparency encoded colors."""
    if weights is not None:
        alphas = weights / weights.max()
        color_r, color_g, color_b = to_rgb(base_color)
        color = [(color_r, color_g, color_b, alpha) for alpha in alphas]
    else:
        color = None

    return color
