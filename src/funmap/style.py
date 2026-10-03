"""
Style sheet and plotting theme for FunMaP publication figures.

Follows APS / Physical Review guidelines:
- Typography: Helvetica / Arial sans-serif
- High-contrast, colorblind-friendly palettes
- Strict bounding boxes (no empty white margins)
- Vector SVG with editable text + 600 DPI raster PNG
"""

import matplotlib as mpl
import matplotlib.pyplot as plt


def apply_publication_style():
    """
    Apply standard publication plot parameters.
    """
    mpl.rcParams.update({
        'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'font.family': 'sans-serif',
        'font.size': 8.5,
        'axes.labelsize': 9.0,
        'axes.titlesize': 9.5,
        'xtick.labelsize': 8.0,
        'ytick.labelsize': 8.0,
        'legend.fontsize': 7.5,
        'figure.titlesize': 10.0,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.5,
        'ytick.minor.width': 0.5,
        'xtick.major.size': 3.5,
        'ytick.major.size': 3.5,
        'xtick.minor.size': 2.0,
        'ytick.minor.size': 2.0,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'lines.linewidth': 1.2,
        'lines.markersize': 5.0,
        'legend.frameon': False,
        'svg.fonttype': 'none',  # Keep text as editable vector text in SVG
        'figure.autolayout': False
    })


# Color definitions
COLORS = {
    'batch1': '#1f77b4',       # Dark blue / primary
    'batch2': '#ff7f0e',       # Amber / orange
    'as_dep': '#7f7f7f',       # Neutral grey
    'single_cap': '#2ca02c',   # Green
    'ideal_sim': '#d62728',    # Red
    'neck': '#1e88e5',         # Blue for metallic necks
    'rim': '#d81b60',          # Red/pink for equatorial rim
    'gold_cap': '#d4af37'      # Metallic gold for caps
}


def save_figure(fig, base_path, dpi=600):
    """
    Save figure in SVG, PDF, and high-resolution PNG with strict tight bounding box.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    base_path : str
        Filepath without extension.
    dpi : int
        Resolution for PNG.
    """
    svg_path = f"{base_path}.svg"
    pdf_path = f"{base_path}.pdf"
    png_path = f"{base_path}.png"
    
    fig.savefig(svg_path, format='svg', bbox_inches='tight', pad_inches=0.04)
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight', pad_inches=0.04)
    fig.savefig(png_path, format='png', dpi=dpi, bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)
