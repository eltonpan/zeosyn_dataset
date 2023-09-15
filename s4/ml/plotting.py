"""Plotting utility functions."""

from typing import Tuple, Dict, Optional, List

import matplotlib.pyplot as plt
import numpy
import pandas
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pandas.api.types import is_float_dtype, is_bool_dtype
from pymatgen.core import Element
from sklearn.metrics import confusion_matrix
from scipy.stats import norm

__all__ = [
    'generate_pt_triangulation',
    'periodic_table_heatmap',
    'periodic_table_heatmap_sbs',
    'scatter_matrix_plot',
]


def skip_7th_row(true_row, inv=False):
    """Utility function to skip the 7th row in periodic table, whose
    chemical elements are rare in the synthesis dataset."""
    if inv:
        if true_row >= 7:
            return true_row + 1
        return true_row
    if true_row >= 7:
        true_row -= 1
    return true_row


def generate_pt_triangulation(
        ncols: int, rows2id: Dict[int, int], padding=(0.1, 0.1, 0.1, 0.1),
        spacing=1.5, header=0.5):
    """
    Generate triangulation for plotting side-by-side periodic tables.

    :param ncols: Number of columns.
    :param rows2id: Periodic-table row id to index-based row id.
    :param padding: Padding of the triangles inside each box.
    :param spacing: Height of each box.
    :param header: Height of the header (i.e., chemical element symbol).
    :returns: Generated triangulation.
    """

    # pylint: disable=too-many-locals
    n_rows_total = len(rows2id)
    y_space0 = numpy.arange(n_rows_total + 1) * spacing - padding[2]
    y_space1 = numpy.arange(n_rows_total + 1) * spacing + header + padding[0]
    x_space0 = numpy.arange(ncols) + padding[3]
    x_space1 = numpy.arange(ncols) + (1 - padding[1])

    x_space = numpy.array([x_space0, x_space1]).T.flatten()
    y_space = numpy.array([y_space0, y_space1]).T.flatten()[:-1]
    x_grid, y_grid = numpy.meshgrid(x_space, y_space)
    x_grid, y_grid = x_grid.flatten(), y_grid.flatten()
    triangles_upper = []
    triangles_lower = []

    n_boxes_a_row = 2 * ncols
    for elem in Element:
        if elem.row not in rows2id:
            continue

        triangles_upper.append([
            n_boxes_a_row * 2 * rows2id[elem.row] + n_boxes_a_row + 2 * (elem.group - 1),
            n_boxes_a_row * 2 * rows2id[elem.row] + n_boxes_a_row + 1 + 2 * (elem.group - 1),
            n_boxes_a_row * 2 * rows2id[elem.row] + n_boxes_a_row * 2 + 2 * (elem.group - 1)
        ])
        triangles_lower.append([
            n_boxes_a_row * 2 * rows2id[elem.row] + n_boxes_a_row + 1 + 2 * (elem.group - 1),
            n_boxes_a_row * 2 * rows2id[elem.row] + n_boxes_a_row * 2 + 2 * (elem.group - 1),
            n_boxes_a_row * 2 * rows2id[elem.row] + n_boxes_a_row * 2 + 2 * (elem.group - 1) + 1
        ])

    return x_grid.flatten(), -y_grid.flatten(), triangles_upper, triangles_lower


def periodic_table_heatmap_sbs(  # pylint: disable=too-many-arguments
        elemental_data1: Dict[str, float], elemental_data2: Dict[str, float],
        cbar_label: str = "",
        cbar_label_size: int = 10,
        show_plot: bool = False,
        cmap: str = "YlOrRd",
        cmap_range: Optional[Tuple[float, float]] = None,
        cax=None,
        blank_color: str = "w",
        value_format: Optional[str] = None,
        include_rows: Optional[List[int]] = None,
        ax=None,
        pt_label_fontsize=9,
        pt_value_fontsize=6,
):
    """
    Plot heatmaps side-by-side.

    :param elemental_data1: The first set of elemental data.
    :param elemental_data2: The second set of elemental data.
    :param cbar_label: Label for the colorbar.
    :param cbar_label_size: Label size for the colorbar.
    :param show_plot: Whether to show this plot after plotting.
    :param cmap: What colormap to use.
    :param cmap_range: Range of the colormap.
    :param cax: The ax to put colorbar.
    :param blank_color: Color of blank elemental data.
    :param value_format: Formatter to use for values.
    :param include_rows: What rows to include in the plot.
    :param ax: Ax to put the periodic table.
    :param pt_label_fontsize: Font size for periodic table cells.
    :param pt_value_fontsize: Font size for periodic table cell values.
    """

    # pylint: disable=too-many-locals
    if cmap_range is None:
        max_val = max(max(elemental_data1.values()), max(elemental_data2.values()))
        min_val = min(min(elemental_data1.values()), min(elemental_data2.values()))
        mean_val = max(abs(max_val), abs(min_val))  # / 2
        min_val, max_val = -mean_val, mean_val
    else:
        min_val, max_val = cmap_range

    if include_rows is None:
        include_rows = {1, 2, 3, 4, 5, 6, 7, 8, 9}
    rows2id = {}
    for i in sorted(include_rows):
        rows2id[i] = len(rows2id)

    value_table = numpy.empty((len(rows2id), 18, 2)) * numpy.nan
    blank_value = min_val - 0.01

    for elem in Element:
        if elem.row not in rows2id:
            continue
        value_table[rows2id[elem.row], elem.group - 1, 0] = elemental_data1.get(
            elem.symbol, blank_value)
        value_table[rows2id[elem.row], elem.group - 1, 1] = elemental_data2.get(
            elem.symbol, blank_value)

    if ax is None:
        ax = plt.gca()

    # We set nan type values to masked values (ie blank spaces)
    data_mask = numpy.ma.masked_invalid(value_table.tolist())

    tri_x, tri_y, tri_upper, tri_lower = generate_pt_triangulation(18, rows2id)

    ax.pcolor(
        numpy.arange(18 + 1),
        numpy.arange(len(rows2id) + 1) * 1.5,
        numpy.full_like(data_mask[:, :, 0], fill_value=float(0)),
        cmap='binary',
        # alpha=0.0,
        edgecolors="k",
        linewidths=1.,
        vmin=0,
        vmax=1,
    )
    heatmap = ax.tripcolor(
        tri_x, -tri_y, tri_upper + tri_lower,
        (
                [data_mask[rows2id[el.row], el.group - 1, 0]
                 for el in Element if el.row in rows2id] +
                [data_mask[rows2id[el.row], el.group - 1, 1]
                 for el in Element if el.row in rows2id]
        ),
        cmap=cmap,
        edgecolors="white",
        linewidths=0.5,
        vmin=min_val - 0.001,
        vmax=max_val + 0.001,
    )

    if cax is None:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='2%', pad=0.1)
    cbar = plt.gcf().colorbar(heatmap, cax=cax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=cbar_label_size)

    # Grey out missing elements in input data
    cbar.cmap.set_under(blank_color)

    # Set the colorbar label and tick marks
    if cbar_label:
        cbar.set_label(cbar_label, rotation=270, labelpad=25, size=cbar_label_size)

    # Refine and make the table look nice
    ax.axis("off")
    ax.invert_yaxis()

    # Label each block with corresponding element and value
    for element in Element:
        if element.row not in rows2id:
            continue

        i, j = rows2id[element.row], element.group - 1
        values = value_table[i, j]

        if numpy.all(numpy.isnan(values)):
            continue

        # print(i, j, symbol)
        ax.text(
            j + 0.5,
            i * 1.5 + 0.4,
            element.symbol,
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=pt_label_fontsize,
            color='k',
        )
        if values[0] != blank_value and value_format is not None:
            ax.text(
                j + 0.12,
                i * 1.5 + 0.62,
                value_format % values[0],
                horizontalalignment="left",
                verticalalignment="center",
                fontsize=pt_value_fontsize,
                color='k',
            )
        if values[1] != blank_value and value_format is not None:
            ax.text(
                j + 0.88,
                i * 1.5 + 1.0,
                value_format % values[1],
                horizontalalignment="right",
                verticalalignment="center",
                fontsize=pt_value_fontsize,
                color='k',
            )

    if show_plot:
        plt.show()


def periodic_table_heatmap(  # pylint:disable=too-many-arguments
        elemental_data: Dict[str, float],
        cbar_label: str = "",
        cbar_label_size: int = 14,
        show_plot: bool = False,
        cmap: str = "YlOrRd",
        cmap_range: Optional[Tuple[float, float]] = None,
        blank_color: str = "w",
        value_format: Optional[str] = None,
        max_row: int = 9,
        ax=None,
):
    """
    Plot heatmap in a periodic table.

    Copied from pymatgen.util.plotting.

    :param elemental_data: The dictionary of elemental data.
    :param cbar_label: Label for the colorbar.
    :param cbar_label_size: Label size for the colorbar.
    :param show_plot: Whether to show this plot after plotting.
    :param cmap: What colormap to use.
    :param cmap_range: Range of the colormap.
    :param cax: The ax to put colorbar.
    :param blank_color: Color of blank elemental data.
    :param value_format: Formatter to use for values.
    :param max_row: Maximal number of rows.
    :param ax: Ax to put the periodic table.
    """

    # pylint:disable=too-many-locals
    max_val, min_val = max(elemental_data.values()), min(elemental_data.values())
    max_row = min(max_row, 9)

    value_table = numpy.empty((max_row, 18)) * numpy.nan
    blank_value = min_val - 0.01

    for el in Element:
        if el.row > max_row:
            continue
        value = elemental_data.get(el.symbol, blank_value)
        value_table[el.row - 1, el.group - 1] = value

    if ax is None:
        ax = plt.gca()

    # We set nan type values to masked values (ie blank spaces)
    data_mask = numpy.ma.masked_invalid(value_table.tolist())
    heatmap = ax.pcolor(
        data_mask,
        cmap=cmap,
        edgecolors="k",
        linewidths=1,
        vmin=min_val - 0.001,
        vmax=max_val + 0.001,
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = plt.gcf().colorbar(heatmap, cax=cax, orientation='vertical')

    #     cbar = ax.colorbar(heatmap)

    # Grey out missing elements in input data
    cbar.cmap.set_under(blank_color)

    # Set the colorbar label and tick marks
    cbar.set_label(cbar_label, rotation=270, labelpad=25, size=cbar_label_size)
    cbar.ax.tick_params(labelsize=cbar_label_size)

    # Refine and make the table look nice
    ax.axis("off")
    ax.invert_yaxis()

    # Label each block with corresponding element and value
    for i, row in enumerate(value_table):
        for j, el in enumerate(row):
            if not numpy.isnan(el):
                symbol = Element.from_row_and_group(i + 1, j + 1).symbol
                ax.text(
                    j + 0.5,
                    i + 0.25,
                    symbol,
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=14,
                    color='k',
                )
                if el != blank_value and value_format is not None:
                    ax.text(
                        j + 0.5,
                        i + 0.5,
                        value_format % el,
                        horizontalalignment="center",
                        verticalalignment="center",
                        fontsize=10,
                        color='k',
                    )

    if show_plot:
        plt.show()


def scatter_matrix_plot(data_frame: pandas.DataFrame,  # pylint:disable=too-many-arguments
                        figsize: Tuple[int, int] = (10, 10),
                        fontsize: int = 9, confusion_matrix_fontsize: int = 9,
                        hist_bins: int = 20,
                        binary_bins: int = 20,
                        scatter_alpha: float = 0.3, scatter_size: int = 2):
    """
    Plot scatter matrix, just as pandas.plotting.scatter_matrix.

    :param data_frame: The data frame to plot.
    :param figsize: Figure size.
    :param fontsize: Fontsize in the plot.
    :param confusion_matrix_fontsize: Fontsize for the labels in the confusion matrix.
    :param hist_bins: Number of bins for histograms.
    :param binary_bins: Number of bins for the binary variable histograms.
    :param scatter_alpha: Alpha of the points in the scatter plot.
    :param scatter_size: Size of the points in the scatter plot.
    """

    # pylint:disable=too-many-branches,too-many-statements,too-many-locals
    ncols = len(data_frame.columns)

    fig, axes = plt.subplots(figsize=figsize, nrows=ncols, ncols=ncols)

    def parse_series(name):
        series = data_frame[name]
        numeric = series.astype(float)
        series_limit = (numeric.min(), numeric.max())
        series_range = series_limit[1] - series_limit[0]
        series_limit = (series_limit[0] - 0.1 * series_range, series_limit[1] + 0.1 * series_range)
        return series, numeric, series_limit, series_range

    def bernoulli_ci(prob, n):
        z1, z2 = norm.ppf(0.025), norm.ppf(0.975)
        return (
            prob - z1 * numpy.sqrt(prob * (1 - prob) / n),
            prob + z1 * numpy.sqrt(prob * (1 - prob) / n)
        )

    def group_by_bins(numeric_series, binary_series):
        binned = pandas.cut(numeric_series, binary_bins)
        probs = data_frame.groupby(binned).apply(
            lambda x: x[binary_series.name].sum() / (len(x) + 1e-5))
        probs_ci = data_frame.groupby(binned).apply(
            lambda x: bernoulli_ci(x[binary_series.name].sum() / (len(x) + 1e-5), len(x) + 1e-5))

        probs.sort_index(inplace=True)
        probs_ci.sort_index(inplace=True)
        bin_locations = sorted([(x.left + x.right) / 2 for x in probs.index])
        return bin_locations, probs, probs_ci

    for i, ax_row in enumerate(axes):
        y_series, y_numeric, y_limit, y_range = parse_series(data_frame.columns[i])

        for j, ax in enumerate(ax_row):
            x_series, x_numeric, x_limit, x_range = parse_series(data_frame.columns[j])

            if i == j:
                ax.hist(x_numeric, bins=hist_bins)
            elif is_bool_dtype(y_series.dtype) and is_float_dtype(x_series.dtype):
                # X is continuous and Y is binary, plot probability
                x_locations, y_probs, y_ci = group_by_bins(x_series, y_series)
                ci_lower, ci_upper = zip(*y_ci.values)
                ax.fill_between(x_locations, ci_lower, ci_upper, color='tab:green', alpha=0.5)
                ax.plot(x_locations, y_probs, 'o-')
            elif is_float_dtype(y_series.dtype) and is_bool_dtype(x_series.dtype):
                # Y is continuous and X is binary, plot probability
                y_locations, x_probs, x_ci = group_by_bins(y_series, x_series)
                ci_lower, ci_upper = zip(*x_ci.values)
                ax.fill_betweenx(y_locations, ci_lower, ci_upper, color='tab:green', alpha=0.5)
                ax.plot(x_probs, y_locations, 'o-')
            elif is_float_dtype(y_series.dtype) and is_float_dtype(x_series.dtype):
                # A normal scatter plot
                ax.scatter(x_series, y_series, s=scatter_size, alpha=scatter_alpha)
            elif is_bool_dtype(x_series.dtype) and is_bool_dtype(y_series.dtype):
                # Both are binary variable, plot confusion matrix
                cmat = confusion_matrix(y_series, x_series)
                ax.imshow(cmat, cmap='summer')
                ax.text(0.2, 0.2, 'TN =\n%d' % cmat[0][0], ha='center', va='center',
                        fontsize=confusion_matrix_fontsize)
                ax.text(0.8, 0.2, 'FP =\n%d' % cmat[0][1], ha='center', va='center',
                        fontsize=confusion_matrix_fontsize)
                ax.text(0.2, 0.8, 'FN =\n%d' % cmat[1][0], ha='center', va='center',
                        fontsize=confusion_matrix_fontsize)
                ax.text(0.8, 0.8, 'TP =\n%d' % cmat[1][1], ha='center', va='center',
                        fontsize=confusion_matrix_fontsize)
                corr = numpy.corrcoef(x_series, y_series)[0][1]
                ax.text(0.5, 0.5, '$\\hat{y}=x$\nCorr = %.3f' % corr,
                        ha='center', va='center', fontsize=confusion_matrix_fontsize)
            else:
                ax.text(sum(x_range) / 2, sum(y_range) / 2,
                        "Don't understand how to plot\n x=%r, y=%r" % (
                            x_series.dtype, y_series.dtype),
                        ha='center', va='center')

            ax.set_xlim(x_limit)
            if i != j:
                ax.set_ylim(y_limit)

            if j > 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel(data_frame.columns[i], fontsize=fontsize)
            if i < len(axes) - 1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(data_frame.columns[j], fontsize=fontsize)

            ax.tick_params(axis='both', which='major', labelsize=fontsize)

    plt.subplots_adjust(wspace=0, hspace=0, left=0.05, right=0.95, top=0.95, bottom=0.05)

    return fig, axes
