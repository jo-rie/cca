from os.path import join, isdir
from os import listdir, makedirs
from typing import List

import matplotlib as mpl
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.patches import Rectangle
import numpy as np
import yaml
from yaml.loader import SafeLoader
from scipy.io import loadmat

from python_src.initialise_plots import update_mpl_rcparams, MPL_ALPHA, MPL_S, plot_style_ca, \
    list_of_cmaps
from python_src.plot_scripts import save_fig, fig_with_size

num_dim = 5  # Number of dimensions
grid_size = 101  # Number of grid points in each dimension

# %%  Plot and directory setup
# Import plot settings
update_mpl_rcparams()
mpl.use('pdf')

# Import yml files
with open('config.yml') as f:
    config = yaml.load(f, Loader=SafeLoader)
with open('config_local.yml') as f:
    config_local = yaml.load(f, Loader=SafeLoader)
turbine_data_root = join(config_local['turbine_base_path'], 'results_copca_20231018')
turbine_plots_root = config_local['turbine_plots_path']


# %% Compute likely empirical profile values for independence copula
# Load Sampled 95%-intervals of profiles for independence copula

def get_empirical_rectangle_from_independence(profile: str, quantile: float = 0.95):
    """
    This function calculates the empirical rectangle from independence for a given profile and quantile.

    Parameters:
    profile (str): The profile to be used. It should be a string that matches the keys in the MATLAB data file.
    quantile (float): The quantile to be used for calculating the rectangle. Default is 0.95.

    Returns:
    Rectangle: A matplotlib.patches.Rectangle object representing the empirical rectangle. The rectangle is defined by its lower left corner, width, and height. The color of the rectangle is grey with an alpha transparency of MPL_ALPHA.

    """
    # Load the MATLAB data file
    matlab_data = loadmat(f'results/bootstrap_independence_profiles/20231026_results_n_{grid_size}.mat')

    # Calculate the lower and upper bounds of the rectangle using the specified quantile
    x_l = np.quantile(matlab_data[f'min_{profile}_profiles'][:, 0], q=1 - quantile)
    y_l = np.quantile(matlab_data[f'min_{profile}_profiles'][:, 1], q=1 - quantile)
    x_r = np.quantile(matlab_data[f'max_{profile}_profiles'][:, 0], q=quantile)
    y_u = np.quantile(matlab_data[f'max_{profile}_profiles'][:, 1], q=quantile)

    # Return a Rectangle object with the calculated bounds and color settings
    return Rectangle((x_l, y_l), x_r - x_l, y_u - y_l, edgecolor=None, facecolor='grey', alpha=MPL_ALPHA)


# Build up list of scenarios
# list_folders = [f for f in listdir(turbine_data_root) if isdir(join(turbine_data_root, f))]
list_folders = ['t_90_p_5_s_230_n_26']

for folder_loop in tqdm(list_folders):
    # print(f'Starting {folder_loop}...', end='')
    path_loop = join(turbine_data_root, folder_loop)
    makedirs(join(turbine_plots_root, folder_loop), exist_ok=True)
    # Compute Axis scaling
    axis_lims = np.zeros(2)
    for i in range(1, num_dim + 1):
        for j in range(1, num_dim + 1):
            if i < j:
                # Load data
                data_matlab = loadmat(join(path_loop, f'dim1_{i}_dim2_{j}.mat'))
                axis_lims[1] = np.nanmax([data_matlab['rowProfiles'][:, [0, 1]].max(),
                                          data_matlab['columnProfiles'][:, [0, 1]].max(),
                                          axis_lims[1]])
                axis_lims[0] = np.nanmin([data_matlab['rowProfiles'][:, [0, 1]].min(),
                                          data_matlab['columnProfiles'][:, [0, 1]].min(),
                                          axis_lims[0]])

    # Plot row and column profiles and data scatter
    for i in range(1, num_dim + 1):
        for j in range(1, num_dim + 1):
            if i < j:
                # Load data
                data_matlab = loadmat(join(path_loop, f'dim1_{i}_dim2_{j}.mat'))

                # Plot row and column profiles
                for (profile, f_or_g) in zip(['row', 'column'], ['F', 'G']):
                    fig, ax = fig_with_size(3, layout='tight', nb_vertically=6)
                    ax.add_artist(get_empirical_rectangle_from_independence(profile))
                    ax.scatter(
                        data_matlab[f'{profile}Profiles'][:, 0].flatten(),
                        data_matlab[f'{profile}Profiles'][:, 1].flatten(),
                        c=np.linspace(0.4, 1, data_matlab[f'{profile}Profiles'].shape[0]),
                        cmap=mpl.colormaps[list_of_cmaps[0]],
                        **plot_style_ca
                    )
                    ax.set(xlabel=f'${f_or_g}_{{:,1}}$', ylabel=f'${f_or_g}_{{:,2}}$', xlim=axis_lims * 1.05, ylim=axis_lims * 1.05)
                    save_fig(fig, join(turbine_plots_root, folder_loop, f'dim1_{i}_dim2_{j}_{profile}Profiles.pdf'))

                # Checkerboard copula plots
                fig, ax = fig_with_size(3, layout='tight', nb_vertically=6)
                matrix = data_matlab['checkerboard_pdf']
                axesim = ax.imshow(
                    matrix,
                    origin='lower', interpolation='nearest', aspect='equal',
                    extent=[0, 1, 0, 1]
                )
                fig.colorbar(axesim)
                save_fig(fig, join(turbine_plots_root, folder_loop, f'dim1_{i}_dim2_{j}_checkerboard.pdf'))

    # print('Finished.')
