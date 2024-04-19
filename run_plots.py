from os.path import join

import matplotlib as mpl
import numpy as np
import yaml
from yaml.loader import SafeLoader

from python_src.initialise_plots import update_mpl_rcparams
from python_src.plot_scripts import parametric_copula_plots, asymmetric_copula_plots, \
    plot_hermite_polys, invalid_copula_plots, comonotonocity_independence_plots, counterexamples_geom_dim_plots, \
    frobenius_distance_plots, fgm_plots, plot_hellinger

# %%  Plot and directory setup
# Import plot settings
update_mpl_rcparams()
mpl.use('pgf')

# Import yml files
with open('config.yml') as f:
    config = yaml.load(f, Loader=SafeLoader)
with open('config_local.yml') as f:
    config_local = yaml.load(f, Loader=SafeLoader)
plot_folder = join(config_local['base_path'], config['plot_dir'])
data_root = join(config_local['base_path'], config['results_dir'])

#%% Run plots
for (copula_short, thetas_to_plot, indices_to_plot, param_name) in zip(
    ['fgm', 'ca', 'amh', 'gumbel', 'gaussian'],
    [
        [.8],  # fgm
        [0.25, 0.5, 0.75],  # ca
        [0.3, 0.5, 0.7, 0.9],  # amh
        [2.5, 5, 7.5, 10],  # gumbel
        [0.25, 0.5, 0.75]  # gaussian
    ],
    [[1]] + 4 * [np.arange(1, 6)],
    [r'\theta', r'\theta', r'\theta', r'\theta', r'\rho']
):
    parametric_copula_plots(
        copula_short=copula_short,
        thetas_to_plot=thetas_to_plot,
        indices_to_plot=indices_to_plot,
        param_name=param_name,
        plot_folder=plot_folder,
        data_root=data_root
    )

fgm_plots(plot_folder=plot_folder, data_root=data_root)
plot_hermite_polys(5, plot_folder=plot_folder, data_root=data_root)

asymmetric_copula_plots(
    indices_to_plot=np.arange(1, 3),
    plot_folder=plot_folder,
    data_root=data_root
)

invalid_copula_plots(
    plot_folder=plot_folder,
    data_root=data_root
)

comonotonocity_independence_plots(
    plot_folder=plot_folder,
    data_root=data_root
)
#
counterexamples_geom_dim_plots(
    plot_folder=plot_folder,
    data_root=data_root
)

frobenius_distance_plots(
    plot_folder=plot_folder,
    data_root=data_root
)

#%% Plots for the appendix
plot_hellinger(
    plot_folder=plot_folder,
    data_root=data_root
)