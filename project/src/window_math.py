"""Window functions from first principles (explicit loops)."""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

WindowName = Literal["rectangular", "hann", "gaussian"]


def rectangular(n: int) -> np.ndarray:
    if n <= 0:
        raise ValueError("Window length must be positive.")
    w = np.empty(n, dtype=float)
    for i in range(n):
        w[i] = 1.0
    return w


def hann(n: int) -> np.ndarray:
    if n <= 0:
        raise ValueError("Window length must be positive.")
    if n == 1:
        return np.array([1.0], dtype=float)
    w = np.empty(n, dtype=float)
    for i in range(n):
        w[i] = 0.5 * (1.0 - math.cos(2.0 * math.pi * i / (n - 1)))
    return w


def gaussian(n: int, sigma: float = 0.4) -> np.ndarray:
    if n <= 0:
        raise ValueError("Window length must be positive.")
    if sigma <= 0 or sigma > 0.5:
        raise ValueError("sigma must be in (0, 0.5].")
    if n == 1:
        return np.array([1.0], dtype=float)
    center = (n - 1) / 2.0
    denom = sigma * center
    w = np.empty(n, dtype=float)
    for i in range(n):
        w[i] = math.exp(-0.5 * ((i - center) / denom) ** 2)
    return w


def get_window(name: WindowName, n: int, *, sigma: float | None = None) -> np.ndarray:
    if name == "rectangular":
        return rectangular(n)
    if name == "hann":
        return hann(n)
    if name == "gaussian":
        if sigma is None:
            return gaussian(n)
        return gaussian(n, sigma=sigma)
    raise ValueError(f"Unsupported window name: {name}")


def window_energy(window: np.ndarray) -> float:
    total = 0.0
    for v in window:
        total += float(v) * float(v)
    return total
