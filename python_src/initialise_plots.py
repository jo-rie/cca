import matplotlib as mpl
import matplotlib.pyplot as plt

MPL_ALPHA = .5
MPL_S = .5

cycler_sym = mpl.cm.get_cmap('PRGn')  # plt.cm.PRGn  # Symmetric plot colors
cycler_01 = mpl.cm.get_cmap('YlGn')  # 0-1 plot colors
plot_style_ca = {'s': MPL_S, 'alpha': MPL_ALPHA, 'vmin': 0, 'vmax': 1}
list_of_cmaps = ['Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']

# Widths and heights
# text_width_pts = 468 # Elsevier
text_width_pts = 441  # De Gruyter
pts_to_inch = 1 / 72.27
text_width = text_width_pts * pts_to_inch
text_height = 661 * pts_to_inch
default_fig_height = text_width / 3.5
default_fig_width = text_width
fig_factor_horizontally = 1.05  # Additional space for room between figures
fig_factor_vertically = 1.4  # Additional space for room between figures for caption etc.


def update_mpl_rcparams():
    """Update matplotlib plot parameters to use pgf and latex."""
    plt.rcParams.update({
        'figure.dpi': 600,
        "text.usetex": True,
        'font.size': 4,
        "font.family": "serif",
        "font.serif": ["Palatino"],
        "figure.figsize": (default_fig_width, default_fig_height),
        'axes.labelsize': 4,
        'legend.fontsize': 4,
    })
    mpl.use('pgf')
