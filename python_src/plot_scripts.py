import os.path
from math import ceil
from os.path import join, exists
from os import makedirs
from typing import Tuple, List, Union

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, patches, ticker
import matplotlib as mpl
from scipy.integrate import quad
from scipy.io import loadmat
from numpy.polynomial.hermite_e import herme2poly
from numpy.polynomial import Polynomial as P
from scipy import stats, linalg
import seaborn as sns

from python_src.initialise_plots import MPL_S, MPL_ALPHA, default_fig_height, default_fig_width, plot_style_ca, \
    list_of_cmaps, fig_factor_horizontally, text_height, fig_factor_vertically

sing_vec_x_label = '$(j - 0.5)/n$'
sing_vec_y_label = '$u_{i, j}$'


def fgm_plots(plot_folder, data_root):
    """Plots for the FGM copula."""
    # Plot of singular vector
    n_array = [100, 1000, 10000]
    theta = .8
    mpl.use('pdf')
    fig_sing_vecs, ax_sing_vecs = fig_with_size(4)
    for n in n_array:
        matlab_data = loadmat(join(data_root, 'fgm', f'raw_n_{n}.mat'))
        theta_i = np.argwhere(np.abs(matlab_data['thetas'].flatten() - theta) < 1e-5).flatten()
        ax_sing_vecs.scatter(matlab_data['plot_grid'][0, :], matlab_data['sing_vecs'][theta_i, :, 0], label=f'$n={n}$',
                             s=MPL_S)
    ax_sing_vecs.set(xlabel=sing_vec_x_label, ylabel='$u_{1, j}$')
    # ax_sing_vecs.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=3) # Create legend below plot with 3 columns
    fig_sing_vecs.legend(loc='outside upper center', ncol=1)
    save_fig(fig_sing_vecs, join(plot_folder, 'fgm', f'fgm_n_comparison_sing_vecs_theta_{theta:.2f}.pdf'))

    # Plot singular value
    theta_array = np.arange(-1, 1.05, 0.1)
    theta_array[10] = 0
    n = 100
    fig, ax = fig_with_size(4)
    matlab_data = loadmat(join(data_root, 'fgm', f'raw_n_{n}.mat'))
    ax.scatter(matlab_data['thetas'].flatten(), matlab_data['sing_vals'][:, 0].flatten(), s=MPL_S)
    ax.set(xlabel=r'$\theta$')
    save_fig(fig, join(plot_folder, 'fgm', f'fgm_sing_vals_n_{n}.pdf'))


def parametric_copula_plots(copula_short: str, thetas_to_plot: list, indices_to_plot: List[Union[List, np.ndarray]],
                            param_name: str, plot_folder: str, data_root: str):
    """Plots for parametric copulas."""
    print(f'Starting {copula_short.upper()}')
    exp_path_tmp = join(plot_folder, copula_short)
    makedirs(exp_path_tmp, exist_ok=True)

    matlab_data = loadmat(join(data_root, copula_short, 'raw_n_50.mat'))

    param_plot_sing_vecs(copula_short, thetas_to_plot, indices_to_plot, matlab_data, exp_path_tmp)
    param_plot_sing_vals(copula_short, thetas_to_plot, indices_to_plot, param_name, matlab_data, exp_path_tmp)
    param_plot_correspondence_analysis(copula_short, thetas_to_plot, indices_to_plot, param_name, matlab_data,
                                       exp_path_tmp, ismar=False)

    # MAR-Plots
    matlab_data = loadmat(join(data_root, copula_short, 'mar_n_50.mat'))

    param_mar_plot_error(copula_short, thetas_to_plot, indices_to_plot, param_name, matlab_data, exp_path_tmp)
    param_mar_plot_eta(copula_short, thetas_to_plot, indices_to_plot, param_name, matlab_data, exp_path_tmp)
    param_plot_correspondence_analysis(copula_short, thetas_to_plot, indices_to_plot, param_name, matlab_data,
                                       exp_path_tmp, ismar=True)

    # Plot geometric dimension
    matlab_data = loadmat(join(data_root, copula_short, "geom_dim_n_theta.mat"))
    thetas = matlab_data['thetas'].flatten()
    ns = matlab_data['ns'].flatten()
    geom_dims = matlab_data['geom_dim']
    fig, ax = fig_with_size(4)
    for (j, theta) in enumerate(thetas):
        ax.scatter(ns, geom_dims[:, j], label=r'$' + param_name + ' = ' + f'{theta}' + r'$', s=MPL_S,
                   alpha=MPL_ALPHA)
    ax.set(xlabel='n', xscale='log', ylabel=r'$\phi_g$', yscale='log')
    ax.legend()
    save_fig(fig, join(plot_folder, copula_short, f'{copula_short}_geom_dim_n_theta.pdf'))


def param_plot_sing_vecs(copula_short: str, thetas_to_plot: list, indices_to_plot: List[Union[List, np.ndarray]],
                         matlab_data: dict, exp_path_tmp: str):
    """Plots of singular vectors for parametric copulas."""
    # Plot singular vectors per theta
    for theta in thetas_to_plot:
        theta_i = np.argwhere(np.abs(matlab_data['thetas'].flatten() - theta) < 1e-5).flatten()
        fig, ax = fig_with_size(4)
        for j in indices_to_plot:
            ax.scatter(
                matlab_data['plot_grid'].flatten(),
                matlab_data['sing_vecs'][theta_i, :, j - 1].flatten(),
                label=r'$u_{' + f'{j}' + r'}$', s=MPL_S, alpha=MPL_ALPHA)
        # ax.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=len(indices_to_plot))
        fig.legend(loc='outside upper center', ncol=3)
        ax.set(xlabel=sing_vec_x_label)
        save_fig(fig, join(exp_path_tmp, f'{copula_short}_param_{theta:.2f}.pdf'))


def param_plot_sing_vals(copula_short: str, thetas_to_plot: list, indices_to_plot: List[Union[List, np.ndarray]],
                         param_name: str, matlab_data: dict, exp_path_tmp: str):
    """Plots of singular values for parametric copulas."""
    # Plot singular values
    fig, ax = fig_with_size(4)
    for theta in thetas_to_plot:
        theta_i = np.argwhere(np.abs(matlab_data['thetas'].flatten() - theta) < 1e-5).flatten()
        ax.scatter(
            np.arange(1, matlab_data['nb_values_to_save'].flatten() + 1),
            matlab_data['sing_vals'][theta_i, :].flatten(),
            label=r'$' + param_name + f'= {theta:.2f}$',
            s=MPL_S, alpha=MPL_ALPHA)
        ax.set(ylabel=r'$s_i$', xlabel=r'$i$', ylim=[0, 1])
    ax.legend()
    save_fig(fig, join(exp_path_tmp, f'{copula_short}_sing_vals.pdf'))


def param_plot_correspondence_analysis(copula_short: str, thetas_to_plot: list,
                                       indices_to_plot: List[Union[List, np.ndarray]],
                                       param_name: str, matlab_data: dict, exp_path_tmp: str,
                                       ismar: bool):
    """Profile pltos for parametric copulas."""
    if ismar:
        mar_short = '_MAR'
    else:
        mar_short = ''
    # Plot Correspondence Analysis
    fig, ax = fig_with_size(4)
    for (theta_enum, theta) in enumerate(thetas_to_plot):
        theta_i = np.argwhere(np.abs(matlab_data['thetas'].flatten() - theta) < 1e-5).flatten()
        ax.scatter(
            matlab_data[f'rowprofiles{mar_short}'][theta_i, :, 0].flatten(),
            matlab_data[f'rowprofiles{mar_short}'][theta_i, :, 1].flatten(),
            label=r'$' + param_name + f'={theta:.1f}' + r'$',
            c=np.linspace(0.4, 1, len(matlab_data[f'rowprofiles{mar_short}'][theta_i, :, 0].flatten())),
            cmap=mpl.colormaps[list_of_cmaps[theta_enum]],
            **plot_style_ca
        )
    # ax.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center')
    fig.legend(loc='outside upper center', ncol=2)
    ax.set(xlabel='$F_{:,1}$ or $G_{:,1}$', ylabel='$F_{:,2}$ or $G_{:,2}$')
    save_fig(fig, join(exp_path_tmp, f'{copula_short}{mar_short}_correspondence_analysis.pdf'))


def param_mar_plot_error(copula_short: str, thetas_to_plot: list, indices_to_plot: List[Union[List, np.ndarray]],
                         param_name: str, matlab_data: dict, exp_path_tmp: str):
    """Error plots for parametric copulas."""
    for approx_order in range(0, 6):
        fig, ax = fig_with_size(4)
        ax.scatter(matlab_data['thetas'].flatten(), matlab_data['frobenius_error_raw'][:, approx_order].flatten(),
                   label=r'raw', s=.5)
        ax.scatter(matlab_data['thetas'].flatten(), matlab_data['frobenius_error_mar'][:, approx_order].flatten(),
                   label=r'MAR', s=.5)
        ax.legend()
        ax.set(xlabel=r'$' + param_name + '$')
        save_fig(fig,
                 join(exp_path_tmp, f'{copula_short}_fro_errors_approx_order_{approx_order}.pdf'))

    # Plot errors over approximation
    for theta in thetas_to_plot:
        theta_i = np.argwhere(np.abs(matlab_data['thetas'].flatten() - theta) < 1e-5).flatten()
        fig, ax = fig_with_size(4)
        ax.scatter(
            matlab_data['approx_orders'].flatten(),
            matlab_data['frobenius_error_raw'][theta_i, :].flatten(),
            label=r'raw', s=MPL_S)
        ax.scatter(
            matlab_data['approx_orders'].flatten(),
            matlab_data['frobenius_error_mar'][theta_i, :].flatten(),
            label=r'MAR', s=MPL_S)
        ax.legend()
        ax.set(xlabel='$n^*$ ')
        save_fig(fig, join(exp_path_tmp, f'{copula_short}_fro_errors_theta_{theta:.2f}.pdf'))


def param_mar_plot_eta(copula_short: str, thetas_to_plot: list, indices_to_plot: List[Union[List, np.ndarray]],
                       param_name: str, matlab_data: dict, exp_path_tmp: str):
    """Eta plots for parametric copulas."""
    for approx_order in range(0, 6):
        fig, ax = fig_with_size(4)
        ax.scatter(matlab_data['thetas'].flatten(), matlab_data['etas'][:, approx_order].flatten(), s=.5)
        # ax.legend()
        ax.set(xlabel=r'$' + param_name + '$', ylabel=f'$\eta$')
        save_fig(fig, join(exp_path_tmp, f'{copula_short}_etas_approx_order_{approx_order}.pdf'))


def asymmetric_copula_plots(indices_to_plot, plot_folder, data_root):
    """Plots for asymmetric copulas."""
    matlab_data = loadmat(join(data_root, 'asymmetric_copula/results_n_100.mat'))

    exp_path_tmp = join(plot_folder, 'asymmetric_copula')
    makedirs(exp_path_tmp, exist_ok=True)

    for i in range(0, len(matlab_data['a_values'].flatten())):
        filename_praefix = f"a_{matlab_data['a_values'].flatten()[i]}_b_{matlab_data['b_values'].flatten()[i]}" \
                           f"_n_{matlab_data['n'].flatten()[0]}"

        # Plot left singular vectors
        fig, ax = fig_with_size(5)
        for j in indices_to_plot:
            ax.scatter(
                matlab_data['plot_grid'].flatten(),
                matlab_data['sing_vecs_left'][i, :, j - 1].flatten(),
                label=r'$u_{' + f'{j}' + '}$', s=MPL_S, alpha=MPL_ALPHA)
        ax.legend()
        ax.set(xlabel=sing_vec_x_label)
        save_fig(fig, join(exp_path_tmp, f'{filename_praefix}_u.pdf'))

        # Plot right singular vectors
        fig, ax = fig_with_size(5)
        for j in indices_to_plot:
            ax.scatter(
                matlab_data['plot_grid'].flatten(),
                matlab_data['sing_vecs_right'][i, :, j - 1].flatten(),
                label=r'$v_{' + f'{j}' + '}$', s=MPL_S, alpha=MPL_ALPHA)
        ax.legend()
        ax.set(xlabel=sing_vec_x_label)
        save_fig(fig, join(exp_path_tmp, f'{filename_praefix}_v.pdf'))

    # Plot singular values
    nb_sing_vals = 5
    fig, ax = fig_with_size(5)
    for i in range(0, len(matlab_data['a_values'].flatten())):
        ax.scatter(
            np.arange(1, nb_sing_vals + 1),
            matlab_data['sing_vals'][i, :nb_sing_vals].flatten(), s=MPL_S,
            label=f"$a = {matlab_data['a_values'].flatten()[i]}$\nb = {matlab_data['a_values'].flatten()[i]}")
    fig.legend(loc='outside upper center')
    ax.set(xlabel='$i$', ylabel='$s_i$')
    save_fig(fig, join(exp_path_tmp, f'sing_vals.pdf'))

    # Plot correspondence analysis
    fig, ax = fig_with_size(4)
    colors = ['C0', 'C1']
    for i in range(0, len(matlab_data['a_values'].flatten())):
        a = matlab_data['a_values'].flatten()[i]
        b = matlab_data['b_values'].flatten()[i]
        col = colors[i]
        F = matlab_data['rowProfiles'][i, :, :]
        G = matlab_data['columnProfiles'][i, :, :]
        ax.scatter(F[:, 0], F[:, 1],
                   c=np.linspace(0.4, 1, len(F[:, 0])), cmap=mpl.colormaps[list_of_cmaps[2 * i]],
                   **plot_style_ca, label=r'$F, a = ' + f'{a:.1f}' + ', b = ' + f'{b:.1f}' + '$')
        ax.scatter(G[:, 0], G[:, 1], marker='s',
                   c=np.linspace(0.4, 1, len(F[:, 0])), cmap=mpl.colormaps[list_of_cmaps[2 * i + 1]],
                   **plot_style_ca, label=r'$G, a = ' + f'{a:.1f}' + ', b = ' + f'{b:.1f}' + '$')
    fig.legend(loc='outside upper center')
    ax.set(xlabel='$F_{:,1}$ or $G_{:,1}$', ylabel='$F_{:,2}$ or $G_{:,2}$')
    save_fig(fig, join(exp_path_tmp, 'asymmetric.pdf'))


def invalid_copula_plots(plot_folder, data_root):
    """Plots of invalid (negative) values."""
    exp_path_tmp = join(plot_folder, 'invalid_values')
    makedirs(exp_path_tmp, exist_ok=True)

    for matlab_data in [loadmat(join(data_root, 'invalid_values', 'Gumbel_theta_2.5_n_50_approx_order_10.mat')),
                        loadmat(join(data_root, 'invalid_values', 'Gumbel_theta_7.5_n_50_approx_order_10.mat'))]:
        matrices_to_plot = [matlab_data['pdf_true'], matlab_data['pdf_approx'],
                            matlab_data['pdf_approx'] - matlab_data['pdf_approx_corrected'],
                            matlab_data['pdf_approx'] < 0]
        matrix_names = ['pdf_true', 'pdf_approx', 'pdf_correction', 'pdf_indicator_invalid']

        overallmax = np.max([matlab_data['pdf_approx'], matlab_data['pdf_true'], matlab_data['pdf_approx_corrected']])
        overallmin = np.min([matlab_data['pdf_approx'], matlab_data['pdf_true'], matlab_data['pdf_approx_corrected']])

        for (matrix, matrix_name) in zip(matrices_to_plot, matrix_names):
            fig, ax = fig_with_size(3, factor_height=.6, layout='tight')
            if matrix_names in ['pdf_true', 'pdf_approx']:
                axesim = ax.imshow(
                    matrix,
                    origin='lower', interpolation='nearest', aspect='equal',
                    extent=[0, 1, 0, 1],
                    vmax=np.ceil(overallmax * 10) / 10, vmin=np.floor(overallmin * 10) / 10)
            else:
                axesim = ax.imshow(matrix, origin='lower', interpolation='nearest', aspect='equal', extent=[0, 1, 0, 1])
            # fig.tight_layout()
            if matrix_names != 'pdf_indicator_invalid':
                fig.colorbar(axesim)
            fig.tight_layout()
            save_fig(fig, join(
                exp_path_tmp,
                f"gumbel_theta_{matlab_data['theta'].flatten()[0]:.1f}_n_{matlab_data['n'].flatten()[0]}_"
                f"approx_order_{matlab_data['k'].flatten()[0]}_{matrix_name}.pdf"))


def comonotonocity_independence_plots(plot_folder, data_root):
    """Plots for the comonotonoicity and independence copulas."""
    exp_path_tmp = join(plot_folder, 'ind_com')
    makedirs(exp_path_tmp, exist_ok=True)
    # Independence
    matlab_data = loadmat(join(data_root, 'ind_com/ind.mat'))
    plot_single_ca(matlab_data['rowprofiles'], join(exp_path_tmp, 'independence.pdf'))

    # Comonotonicity
    matlab_data = loadmat(join(data_root, 'ind_com/com_mar.mat'))
    plot_single_ca(matlab_data['rowprofiles'], join(exp_path_tmp, 'comonotonicity.pdf'))
    plot_single_ca(matlab_data['rowprofiles_MAR'],
                   join(exp_path_tmp, 'comonotonicity_MAR.pdf'))


def counterexamples_geom_dim_plots(plot_folder, data_root):
    """Plots of the counterexample for the geometric dimension."""
    exp_path_tmp = join(plot_folder, 'counterexample_geometric_dimension')
    makedirs(exp_path_tmp, exist_ok=True)
    plot_rectangles(
        lowerleft_points=[[0.000, 0.000], [0.125, 0.125], [0.500, 0.125], [0.625, 0.000], [0.375, 0.250],
                          [0.750, 0.250],
                          [0.250, 0.375], [0.875, 0.375], [0.125, 0.500], [0.500, 0.500], [0.000, 0.625],
                          [0.625, 0.625],
                          [0.250, 0.750], [0.875, 0.750], [0.375, 0.875], [0.750, 0.875]],
        rec_width=1 / 8,
        filename=join(exp_path_tmp, 'pdf_copula.pdf')
    )
    plot_rectangles(
        lowerleft_points=[[0, 0], [.5, 0],
                          [1 / 4, 1 / 4], [3 / 4, 1 / 4],
                          [0, 1 / 2], [1 / 2, 1 / 2],
                          [1 / 4, 3 / 4], [3 / 4, 3 / 4]],
        rec_width=1 / 4,
        filename=join(exp_path_tmp, 'pdf_disc_copula_n4.pdf')
    )
    # Plot Geometric Dimension
    matlab_data = loadmat(join(data_root, 'counterexample_geometric_dimension', 'geom_dims.mat'))

    fig, ax = fig_with_size(3, factor_height=.75)
    ax.plot(matlab_data['ns'].flatten()[1:],
            matlab_data['geom_dims'].flatten()[1:])

    ax.xaxis.set_major_locator(ticker.MultipleLocator(5.00))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1.00))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())

    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.00))
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

    ax.set(xlim=[0, 8], ylim=[0, 4], xlabel=r'$n$', ylabel=r'$\gamma(\mathbf{ C}^n)$')
    # fig.tight_layout()
    # fig.savefig(join(exp_path_tmp, 'geometric_dimension_over_n.pdf'))
    save_fig(fig, join(exp_path_tmp, 'geometric_dimension_over_n.pdf'))


def frobenius_distance_plots(plot_folder, data_root):
    """Plots of the Frobenius distance as similarity measure."""
    # Plot Frobenius distance for Gaussian copulas with different rhos
    exp_path_tmp = join(plot_folder, 'frobenius_distance')
    makedirs(exp_path_tmp, exist_ok=True)
    matlab_data = loadmat(join(data_root, 'frobenius_distance', 'distance_n_100_rhos.mat'))
    vmin = 0
    vmax = None
    for mat_name in ['results_fro', 'results_d2', 'results_d1']:
        fig, ax = fig_with_size(4, factor_height=.6)
        data = pd.DataFrame(
            matlab_data[mat_name],
            columns=[f'{r:.3f}' for r in matlab_data['rhos'].flatten().tolist()],
            index=[f'{r:.3f}' for r in matlab_data['rhos'].flatten().tolist()])
        sns.heatmap(data, ax=ax, vmin=vmin, vmax=vmax)
        ax.set(xlabel=r'$\rho$', ylabel=r'$\rho$')
        save_fig(fig, join(exp_path_tmp,
                           f'Gaussian_frobenius_dist_n_100_rhos_{mat_name}.pdf'))

    #% Plot Frobenius distance for Gaussian copulas with increasing truncation
    matlab_data = loadmat(join(data_root, 'frobenius_distance', 'Gaussian_over_k.mat'))

    tau1 = matlab_data['tau1'].flatten()
    tau2 = matlab_data['tau2'].flatten()
    k_max_plot = 25
    x_range = np.arange(1, k_max_plot + 1)
    val_d1 = matlab_data['val_d1']
    val_d2 = matlab_data['val_d2']

    fig, ax = fig_with_size(2, factor_height=.7)
    for i in range(len(tau1)):
        ax.plot(x_range, val_d1[i, :k_max_plot], '.',
                linestyle='dashed', markersize=2, linewidth=.5, alpha=.5,
                label=r'$\delta_1, ' + f'{tau1[i]:.1f}' + r', ' + f'{tau2[i]:.1f}$')
        ax.plot(x_range, val_d2[i, :k_max_plot], '.',
                linestyle='dashed', markersize=2, linewidth=.5, alpha=.5,
                label=r'$\delta_2, ' + f'{tau1[i]:.1f}' + r', ' + f'{tau2[i]:.1f}$')
    ax.set(xlabel='$k^\star$')
    # ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left')
    fig.legend(loc='outside upper center', ncol=4, fontsize=4)
    save_fig(fig, join(exp_path_tmp, 'Gaussian_over_k.pdf'))

    #% Plot frobenius distance for different families
    matlab_data = loadmat(join(data_root, 'frobenius_distance', 'frobenius_different_families_n_100.mat'))
    vmin = 0
    vmax = None

    fam_tau = ["G, 0.3", "AMH, 0.3", "Cl, 0.3", "Ga, 0.3", "G, 0.9", "Cl, 0.9", "Ga, 0.9"]

    for mat_name in ['results_fro', 'results_d1', 'results_d2']:
        fig, ax = fig_with_size(4, factor_height=.6)
        data = pd.DataFrame(matlab_data[mat_name], columns=fam_tau,
                            index=fam_tau)
        sns.heatmap(data[[c for c in data.columns if not c.startswith('AMH')]], ax=ax, vmin=vmin, vmax=vmax)
        # sns.heatmap(data, ax=ax, vmin=vmin, vmax=vmax)
        save_fig(fig, join(exp_path_tmp, f'copula_families_{mat_name}.pdf'))

    #% Plot frobenius distance for two copula combinations with increasing grid size
    matlab_data = loadmat(join(data_root, 'frobenius_distance', 'frobenius_increasing_n.mat'))

    cop_str_array = ["Gaussian0.1;Gaussian0.9", "Gaussian0.4,Gaussian0.6", "Gumbel1;Gumbel10", "Gumbel4;Gumbel6"]

    for (i_cop, cop) in enumerate(cop_str_array):
        fig, ax = fig_with_size(4, factor_height=.6)
        n_array = matlab_data['n_array'].flatten()
        ax.loglog(n_array, matlab_data['d1_array'][i_cop, :].flatten(), label=r'$\delta_1$')
        ax.loglog(n_array, matlab_data['d2_array'][i_cop, :].flatten(), label=r'$\delta_2$')
        ax.set(xlabel='$n$')
        # ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left')
        fig.legend(loc='outside upper center')
        save_fig(fig, join(exp_path_tmp, f'increasing_n_{cop}.pdf'))


def plot_hellinger(plot_folder, data_root):
    """Plots for the Hellinger distance."""
    exp_path_tmp = join(plot_folder, 'hellinger')
    makedirs(exp_path_tmp, exist_ok=True)
    matlab_data = loadmat(join(data_root, 'hellinger', 'Gaussian_theta_0.75_n_100.mat'))

    # Plot singular vectors
    for postfix in ['fro', 'hellinger']:
        fig, ax = fig_with_size(2)
        for j in range(5):
            ax.scatter(
                matlab_data['plot_grid'].flatten(),
                matlab_data[f'sing_vecs_{postfix}'][:, j - 1].flatten(),
                label=r'$u_{' + f'{j}' + r'}$', s=MPL_S, alpha=MPL_ALPHA)
        # ax.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=len(indices_to_plot))
        fig.legend(loc='outside upper center', ncol=3)
        ax.set(xlabel=sing_vec_x_label)
        save_fig(fig, join(exp_path_tmp, f'sing_vec_{postfix}.pdf'))

    # Plot singular values
    fig, ax = fig_with_size(3)
    ax.scatter(
        np.arange(1, 6),
        matlab_data['sing_vals_fro'].flatten(),
        label='Frobenius', s=MPL_S)
    ax.scatter(
        np.arange(1, 6),
        matlab_data['sing_vals_hellinger'].flatten(),
        label='Hellinger', s=MPL_S)
    fig.legend()
    ax.set(xlabel='$i$', ylabel='$s_i$')
    save_fig(fig, join(exp_path_tmp, f'sing_vals.pdf'))

    # Plot checkerboards
    for postfix in ['fro', 'hellinger']:
        fig, ax = fig_with_size(3)
        axesim = ax.imshow(
            matlab_data[f'checkerboard_{postfix}'],
            origin='lower', interpolation='nearest', aspect='equal',
            extent=[0, 1, 0, 1]
        )
        fig.colorbar(axesim)
        save_fig(fig, join(exp_path_tmp, f'checkerboard_{postfix}.pdf'))


def plot_rectangles(lowerleft_points, rec_width, filename):
    """Plots of rectangles for the profile plots."""
    fig, ax = fig_with_size(3, factor_height=.75)
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_locator(ticker.MultipleLocator(0.5))
        axis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.set(xlim=[0, 1], ylim=[0, 1], xlabel='u', ylabel='v')
    ax.grid(True, linestyle='--', linewidth=0.2)
    for point in lowerleft_points:
        ax.add_patch(patches.Rectangle(point, width=rec_width, height=rec_width))
    save_fig(fig, filename)


def plot_hermite_polys(max_order, plot_folder, data_root):
    """Plots Hermite polynomials."""
    x_grid = np.linspace(0, 1, 100, endpoint=False) + 1 / 200
    fig, ax = fig_with_size(4)
    for degree in range(1, max_order + 1):
        poly = get_transformed_hermite_poly(degree)
        y = poly(x_grid)
        y = y / linalg.norm(y)
        ax.plot(x_grid, y, label=r'$\psi_' + f'{degree}' + r'(x)$', linewidth=.5)
    ax.set(xlabel='x')
    # ax.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=ceil(max_order / 2))
    fig.legend(loc='outside upper center', ncol=2)
    save_fig(fig, join(plot_folder, 'gaussian', f'transformed_hermite_polys.pdf'))

    degree_pairs = [(0, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4)]
    ns = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]

    prod_by_section_filename = join(data_root, 'gaussian', 'data_hermite_scalar_prod_by_section.npz')
    if os.path.exists(prod_by_section_filename):
        npz_array = np.load(prod_by_section_filename)
        values = npz_array['values']
    else:
        values = np.zeros([len(degree_pairs), len(ns)])
        # Compute data scalar product of data over n
        for (i_d, d) in enumerate(degree_pairs):
            for (i_n, n) in enumerate(ns):
                values[i_d, i_n] = rangewise_scalar(d[0], d[1], n)
        np.savez(prod_by_section_filename, ns=ns, degree_pairs=degree_pairs,
                 values=values)

    fig, ax = fig_with_size(4)
    for (i_d, d) in enumerate(degree_pairs):
        ax.plot(ns[2:], values[i_d, 2:], label=d, marker='o')
    ax.legend()
    save_fig(fig, join(plot_folder, 'gaussian', 'hermite_scalarprod_over_n.pdf'))

    # Compute data of distance between singular vectors and singular functions
    ns = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
    rhos = [.25, .5, .75]
    degrees = np.arange(1, 8, step=1, dtype=int)
    diffs = np.zeros([len(rhos), len(degrees), len(ns)])

    def filename_rangewise_integral(rho, n, degree):
        return join(data_root, 'gaussian_rangewise_integral', f'rho_{rho:.2f}_n_{n}_degree_{d}.npy')

    makedirs(join(data_root, 'gaussian_rangewise_integral'), exist_ok=True)
    l2_dist_filename = join(data_root, 'l2_dist_singular_vectors_hermite_polys.npz')
    if os.path.exists(l2_dist_filename):
        npz_array = np.load(l2_dist_filename)
        diffs = npz_array['diffs']
    else:
        for (i_rho, rho) in enumerate(rhos):
            # tau = taus[1]
            for (i_n, n) in enumerate(ns):
                sing_vecs = loadmat(join(data_root, 'gaussian', f'gaussian_sing_vecs_for_n_{n}_rho_{rho:.2f}.mat'))[
                    'singular_vectors']
                for (i_d, d) in enumerate(degrees):
                    # Compute vector
                    if not (exists(filename_rangewise_integral(rho, n, d))):
                        vec_herme = rangewise_integral_herme(d, n)
                        vec_herme = vec_herme / np.linalg.norm(vec_herme) * np.sign(vec_herme[0])  # normalize
                        np.save(filename_rangewise_integral(rho, n, d), vec_herme)
                    else:
                        vec_herme = np.load(filename_rangewise_integral(rho, n, d))
                    # Compute difference
                    diffs[i_rho, i_d, i_n] = np.min(
                        [np.linalg.norm(vec_herme - sing_vecs[:, d - 1]),
                         np.linalg.norm(vec_herme + sing_vecs[:, d - 1])])
        np.savez(l2_dist_filename, ns=ns, degrees=degrees, taus=rhos, diffs=diffs)

    for (i_rho, rho) in enumerate(rhos):
        fig, ax = fig_with_size(3, factor_height=.6)
        for (i_d, d) in enumerate(degrees):
            ax.loglog(ns, diffs[i_rho, i_d, :], marker='o', linestyle='dashed', label=f'd={d}', markersize=MPL_S,
                      linewidth=.2)
        # ax.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=len(degrees))
        fig.legend(loc='outside upper center', ncol=4)
        ax.set(xlabel='$n$')
        save_fig(fig, join(plot_folder, 'gaussian',
                           f'gaussian_decomposition_distance_hermite_integration_by_parts_rho_{rho:.2f}.pdf'))


def fig_with_size(nb_horizontally=1, nb_vertically=1, fig_height=None,
                  fig_width=None, factor_height=None, layout='constrained') -> Tuple[plt.Figure, plt.Axes]:
    """Return a figure so that nb_horizontally fit next to each other and nb_horizontally fit below each other"""
    if fig_height is None:
        if factor_height is None:
            if nb_vertically == 1:
                fig_height = default_fig_height
            else:
                fig_height = text_height / (nb_vertically * fig_factor_vertically)
        else:
            fig_height = default_fig_height * factor_height
    if fig_width is None:
        fig_width = default_fig_width / (nb_horizontally * fig_factor_horizontally)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), layout=layout)
    return fig, ax


def save_fig(fig: plt.Figure, path: str):
    """Helper function to save and close a figure."""
    # fig.tight_layout() # Use constrained layout
    # fig.savefig(path, bbox_inches=('tight'))
    fig.savefig(path)
    plt.close(fig)


def get_poly_from_herme(degree: int):
    """Generate a Hermite polynomial of a given degree."""
    return P(herme2poly(degree * [0] + [1]))


def get_transformed_hermite_poly(degree: int):
    """Get a (0,1)-transformed Hermite polynomial."""
    poly = get_poly_from_herme(degree)

    def poly_t(x) -> np.ndarray: return poly(stats.norm.ppf(x))

    return poly_t


def plot_single_ca(mat, filename):
    """Plot a single correspondence analysis plot."""
    fig, ax = fig_with_size(4)
    ax.scatter(
        mat[:, 0], mat[:, 1],
        c=np.linspace(0.4, 1, len(mat[:, 0])),
        cmap=mpl.colormaps[list_of_cmaps[0]],
        **plot_style_ca)
    adjust_lims(ax, 'x')
    adjust_lims(ax, 'y')
    ax.set(xlabel='$F_{:,1}$ or $G_{:,1}$', ylabel='$F_{:,2}$ or $G_{:,2}$')
    ax.locator_params(axis='both', nbins=2)
    save_fig(fig, filename)


def adjust_lims(ax, axis='x', lower_bound=-.01, upper_bound=.01):
    """Helper function to adjust axis limits."""
    if axis == 'x':
        current_lim_lower, current_lim_upper = ax.get_xlim()
    elif axis == 'y':
        current_lim_lower, current_lim_upper = ax.get_ylim()
    else:
        return False
    lim_lower = np.min([current_lim_lower, lower_bound])
    lim_upper = np.max([current_lim_upper, upper_bound])
    if axis == 'x':
        ax.set_xlim(lim_lower, lim_upper)
    if axis == 'y':
        ax.set_ylim(lim_lower, lim_upper)
    return ax


def herme_poly(degree: int):
    """Get a Hermite polynomial of a single degree."""
    return P(herme2poly(degree * [0] + [1]))


def trafo_herme_poly(degree: int):
    """Get a signle transformed Hermite polynomial."""
    poly = herme_poly(degree)

    def poly_t(x) -> np.ndarray: return poly(stats.norm.ppf(x))

    return poly_t


def rangewise_integral_herme(degree: int, n: int) -> np.ndarray:
    """Compute the rangewise integral for a Hermite polynomial."""
    poly_t = trafo_herme_poly(degree)

    av_values = np.zeros(n)

    for i in range(0, len(av_values)):
        # Compute value for poly1
        av_values[i] = quad(lambda x: poly_t(x), i / n, (i + 1) / n)[0]

    return av_values * n


def rangewise_scalar(degree1: int, degree2: int, n: int) -> np.ndarray:
    """Compute the rangewise integral of the scalar product of two Hermite polynomials."""
    av_values1 = rangewise_integral_herme(degree1, n)
    av_values2 = rangewise_integral_herme(degree2, n)
    return ((av_values1 / (np.linalg.norm(av_values1))) * (av_values2 / (np.linalg.norm(av_values1)))).sum()
