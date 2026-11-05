"""ヒルベルト曲線による点の並べ替え。

セルを母点のヒルベルト曲線上の位置（ヒルベルト指数）の順に並べると、番号の近いセルが
空間的にも近くなる。セル単位のループでメモリアクセスが局所的になり、番号の連続した範囲で
領域を分けると、まとまった形の部分領域になる（将来の並列化に使う）。

指数の計算は、2^p × 2^p の格子上の整数座標 (x, y) から、上位のビットから順に象限を
決めて回転・反転する標準的な方法（Hilbert, 1891 の曲線の離散版）による。
"""

from __future__ import annotations

import numpy as np

DEFAULT_ORDER = 16
"""ヒルベルト曲線の次数 p の既定値。領域の外接正方形を 2^p × 2^p に分ける。"""


def hilbert_index(ix: np.ndarray, iy: np.ndarray, order: int = DEFAULT_ORDER) -> np.ndarray:
    """2^order × 2^order の格子上の整数座標 (ix, iy) のヒルベルト指数を返す。

    指数は 0 から 4^order − 1 までの整数で、(0, 0) が 0、(2^order − 1, 0) が最後になる。
    """
    if not 1 <= order <= 31:
        raise ValueError(f"order must be in 1..31, got {order}")
    n = 1 << order
    x = np.asarray(ix, dtype=np.int64).copy()
    y = np.asarray(iy, dtype=np.int64).copy()
    if x.shape != y.shape:
        raise ValueError("ix and iy must have the same shape")
    if x.size and (x.min() < 0 or y.min() < 0 or x.max() >= n or y.max() >= n):
        raise ValueError(f"coordinates must be in 0..{n - 1}")

    d = np.zeros_like(x)
    s = n >> 1
    while s > 0:
        rx = (x & s) > 0
        ry = (y & s) > 0
        d += s * s * ((3 * rx.astype(np.int64)) ^ ry.astype(np.int64))
        # 象限に応じて回転・反転し、次のビットで同じ向きの曲線になるようにする
        flip = ~ry & rx
        x = np.where(flip, n - 1 - x, x)
        y = np.where(flip, n - 1 - y, y)
        swap = ~ry
        x, y = np.where(swap, y, x), np.where(swap, x, y)
        s >>= 1
    return d


def hilbert_order(points: np.ndarray, order: int = DEFAULT_ORDER) -> np.ndarray:
    """点をヒルベルト曲線の順に並べる添字（`points[perm]` が並べ替えた点）を返す。

    点の外接正方形（x と y で同じ縮尺）を 2^order × 2^order に分けて指数を計算する。
    縮尺を揃えるので、細長い領域でも曲線はゆがまない。指数が同じ点は元の順を保つ。
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError("points must have shape (n, 2)")
    if len(pts) == 0:
        return np.zeros(0, dtype=np.int64)
    lower = pts.min(axis=0)
    extent = float((pts.max(axis=0) - lower).max())
    n = 1 << order
    scale = (n - 1) / extent if extent > 0.0 else 0.0
    grid = np.clip(np.floor((pts - lower) * scale), 0, n - 1).astype(np.int64)
    index = hilbert_index(grid[:, 0], grid[:, 1], order)
    return np.argsort(index, kind="stable")
