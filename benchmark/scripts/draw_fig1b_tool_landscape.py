#!/usr/bin/env python3
"""
Fig 1B - Tool Landscape: Technical Approaches Comparison

Usage (in Docker):
    python3 /benchmark/scripts/draw_fig1b_tool_landscape.py
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from pathlib import Path
import os, sys

print(f"Python: {sys.version}")
print(f"matplotlib: {matplotlib.__version__}")
print(f"CWD: {os.getcwd()}")

# Paths - try Docker path first, then local
if Path('/benchmark/icons').exists():
    ICONS_DIR = Path('/benchmark/icons')
    FIGURES_DIR = Path('/benchmark/figures/main')
else:
    ICONS_DIR = Path('benchmark/icons')
    FIGURES_DIR = Path('benchmark/figures/main')

FIGURES_DIR.mkdir(parents=True, exist_ok=True)
print(f"ICONS_DIR: {ICONS_DIR} (exists={ICONS_DIR.exists()})")
print(f"FIGURES_DIR: {FIGURES_DIR} (exists={FIGURES_DIR.exists()})")
print(f"Icons found: {list(ICONS_DIR.glob('*.png'))}")

# Colors
CC_PINK = '#E91E63'
CC_BG = '#FCE4EC'
CARD_BG = '#F5F5F5'
CARD_BORDER = '#BDBDBD'
TXT1 = '#212121'
TXT2 = '#616161'
BADGE_BG = '#FAFAFA'
BADGE_BD = '#9E9E9E'

TOOLS = [
    {
        'name': 'Zellkonverter',
        'desc': 'R\u2192Python bridge (basilisk) \u2014 needs R + Python',
        'icons': ['r_logo.png', 'python_logo.png'],
        'direction': '\u2194',
        'deps': 'Runtime Deps: R, Python',
        'highlight': False,
    },
    {
        'name': 'anndataR',
        'desc': 'Native R (hdf5r) \u2014 needs R only',
        'icons': ['r_logo.png', 'file_logo.png'],
        'direction': '\u2190',
        'deps': 'Runtime Deps: R',
        'highlight': False,
    },
    {
        'name': 'convert2anndata',
        'desc': 'R\u2192Python bridge (reticulate) \u2014 needs R + Python',
        'icons': ['r_logo.png', 'python_logo.png'],
        'direction': '\u2192',
        'deps': 'Runtime Deps: R, Python',
        'highlight': False,
    },
    {
        'name': 'easySCF',
        'desc': 'Custom HDF5 intermediate \u2014 needs R + Python',
        'icons': ['hdf5_logo.png', 'file_logo.png'],
        'direction': '\u2194',
        'deps': 'Runtime Deps: R, Python',
        'highlight': False,
    },
    {
        'name': 'CrossCell',
        'desc': 'Native Rust (no runtime) \u2014 standalone binary',
        'icons': ['rust_logo.png', 'file_logo.png'],
        'direction': '\u2194',
        'deps': 'Runtime Deps: None (Standalone)',
        'highlight': True,
    },
]


def load_icon(filename, target_h=14):
    """Load icon, normalize to same height."""
    path = ICONS_DIR / filename
    if not path.exists():
        return None
    try:
        img = mpimg.imread(str(path))
        # 根据原始高度计算 zoom，使所有图标显示高度一致
        h = img.shape[0]
        zoom = target_h / h
        return OffsetImage(img, zoom=zoom)
    except Exception as e:
        print(f"  ERROR: {e}")
        return None



def draw():
    W = 10.0
    ROW_H = 0.85
    N = len(TOOLS)
    TITLE_H = 0.55
    FIG_H = TITLE_H + N * ROW_H + 0.45

    fig, ax = plt.subplots(figsize=(W, FIG_H))
    ax.set_xlim(0, W)
    ax.set_ylim(0, FIG_H)
    ax.axis('off')
    fig.patch.set_facecolor('white')

    # Title
    ax.text(W / 2, FIG_H - 0.18,
            'Tool Comparison Overview: Technical Approaches',
            fontsize=14, fontweight='bold', ha='center', va='center',
            color=TXT1, fontfamily='sans-serif')

    # Outer frame
    pad = 0.2
    outer = mpatches.FancyBboxPatch(
        (pad, 0.15), W - 2 * pad, FIG_H - TITLE_H - 0.05,
        boxstyle='round,pad=0.1',
        facecolor='white', edgecolor=CARD_BORDER, linewidth=1.2)
    ax.add_patch(outer)

    CL = 0.4
    CW = W - 0.8
    y_start = FIG_H - TITLE_H - 0.12

    for i, t in enumerate(TOOLS):
        yt = y_start - i * ROW_H
        yc = yt - ROW_H / 2
        cy = yt - ROW_H + 0.08

        if t['highlight']:
            bg, bc, bw = CC_BG, CC_PINK, 2.0
        else:
            bg, bc, bw = CARD_BG, CARD_BORDER, 1.0

        card = mpatches.FancyBboxPatch(
            (CL, cy), CW, ROW_H - 0.16,
            boxstyle='round,pad=0.08',
            facecolor=bg, edgecolor=bc, linewidth=bw)
        ax.add_patch(card)

        # Tool name
        nx = CL + 0.25
        ax.text(nx, yc + 0.13, t['name'],
                fontsize=11.5, fontweight='bold', ha='left', va='center',
                color=TXT1, fontfamily='sans-serif')

        # Icons
        for j, ic_file in enumerate(t['icons']):
            im = load_icon(ic_file )
            if im is not None:
                ab = AnnotationBbox(im, (nx + j * 0.38, yc - 0.11),
                                    frameon=False, box_alignment=(0, 0.5))
                ax.add_artist(ab)

        # Description
        dx = nx + len(t['icons']) * 0.38 + 0.18
        ax.text(dx, yc - 0.11, t['desc'],
                fontsize=8.5, ha='left', va='center',
                color=TXT2, fontfamily='sans-serif')

        # Direction badge
        ax_r = CL + CW - 3.2
        db = mpatches.FancyBboxPatch(
            (ax_r - 0.18, yc - 0.14), 0.36, 0.28,
            boxstyle='round,pad=0.05',
            facecolor=BADGE_BG, edgecolor=BADGE_BD, linewidth=0.8)
        ax.add_patch(db)
        ax.text(ax_r, yc, t['direction'],
                fontsize=14, ha='center', va='center',
                color=TXT1, fontfamily='sans-serif')

        # Runtime deps badge
        rx = CL + CW - 1.5
        dtxt = t['deps']
        bw2 = max(len(dtxt) * 0.082, 2.1)
        rb = mpatches.FancyBboxPatch(
            (rx - bw2 / 2, yc - 0.13), bw2, 0.26,
            boxstyle='round,pad=0.05',
            facecolor=BADGE_BG, edgecolor=BADGE_BD, linewidth=0.8)
        ax.add_patch(rb)
        ax.text(rx, yc, dtxt,
                fontsize=8, ha='center', va='center',
                color=TXT2, fontfamily='sans-serif')

    # Save
    for fmt in ['pdf', 'png']:
        out = FIGURES_DIR / f'fig1b_tool_landscape.{fmt}'
        print(f"Saving {out} ...")
        fig.savefig(str(out), dpi=300, bbox_inches='tight', facecolor='white')
        if out.exists():
            print(f"  OK: {out} ({out.stat().st_size} bytes)")
        else:
            print(f"  FAILED: {out} does not exist!")
    plt.close(fig)
    print('Done!')


if __name__ == '__main__':
    draw()
