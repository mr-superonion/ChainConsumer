# -*- coding: utf-8 -*-
import logging
import numpy as np
from scipy.integrate import simps
from scipy.signal import fftconvolve
from scipy.interpolate import interp1d

from .helpers import (
    get_grid_bins,
    get_latex_table_frame,
    std_weight,
    get_extents,
)


class Analysis(object):

    summaries = [
        "max",
        "mean",
        "cumulative",
        "max_symmetric",
        "max_shortest",
        "max_central",
    ]

    def __init__(self, parent):
        self.parent = parent
        self._logger = logging.getLogger("chainconsumer")

        self._summaries = {
            "max": self.get_parameter_summary_max,
            "mean": self.get_parameter_summary_mean,
            "cumulative": self.get_parameter_summary_cumulative,
            "max_symmetric": self.get_paramater_summary_max_symmetric,
            "max_shortest": self.get_parameter_summary_max_shortest,
            "max_central": self.get_parameter_summary_max_central,
        }

    def get_latex_table(
        self,
        parameters=None,
        transpose=False,
        caption=None,
        label="tab:model_params",
        hlines=True,
        blank_fill="--",
        filename=None,
    ):  # pragma: no cover
        """Generates a LaTeX table from parameter summaries.

        Parameters
        ----------
        parameters : list[str], int optional
            A list of what parameters to include in the table. By default, includes all parameters
        transpose : bool, optional
            Defaults to False, which gives each column as a parameter, each chain (framework)
            as a row. You can swap it so that you have a parameter each row and a framework
            each column by setting this to True
        caption : str, optional
            If you want to generate a caption for the table through Python, use this.
            Defaults to an empty string
        label : str, optional
            If you want to generate a label for the table through Python, use this.
            Defaults to an empty string
        hlines : bool, optional
            Inserts ``\\hline`` before and after the header, and at the end of table.
        blank_fill : str, optional
            If a framework does not have a particular parameter, will fill that cell of
            the table with this string.
        filename : str, optional
            The file to save the output string to

        Returns
        -------
        str
            the LaTeX table.
        """
        if parameters is None:
            parameters = self.parent._all_parameters
        elif isinstance(parameters, int):
            parameters = self.parent._all_parameters[:parameters]
        for p in parameters:
            assert isinstance(
                p, str
            ), "Generating a LaTeX table requires all parameters have labels"
        num_parameters = len(parameters)
        chains = self.parent.get_mcmc_chains()
        num_chains = len(chains)
        fit_values = self.get_summary(squeeze=False, chains=chains)
        if label is None:
            label = ""
        if caption is None:
            caption = ""

        end_text = " \\\\ \n"
        if transpose:
            column_text = "c" * (num_chains + 1)
        else:
            column_text = "c" * (num_parameters + 1)

        center_text = ""
        hline_text = "\\hline\n"
        if hlines:
            center_text += hline_text + "\t\t"
        if transpose:
            center_text += (
                " & ".join(["Parameter"] + [c.name for c in chains]) + end_text
            )
            if hlines:
                center_text += "\t\t" + hline_text
            for p in parameters:
                arr = ["\t\t" + p]
                for chain_res in fit_values:
                    if p in chain_res:
                        arr.append(self.get_parameter_text(*chain_res[p], wrap=True))
                    else:
                        arr.append(blank_fill)
                center_text += " & ".join(arr) + end_text
        else:
            center_text += " & ".join(["Model"] + parameters) + end_text
            if hlines:
                center_text += "\t\t" + hline_text
            for name, chain_res in zip([c.name for c in chains], fit_values):
                arr = ["\t\t" + name]
                for p in parameters:
                    if p in chain_res:
                        arr.append(self.get_parameter_text(*chain_res[p], wrap=True))
                    else:
                        arr.append(blank_fill)
                center_text += " & ".join(arr) + end_text
        if hlines:
            center_text += "\t\t" + hline_text
        final_text = get_latex_table_frame(caption, label) % (column_text, center_text)

        if filename is not None:
            with open(filename, "w") as f:
                f.write(final_text)

        return final_text

    def get_summary(self, squeeze=True, parameters=None, chains=None):
        """Gets a summary of the marginalised parameter distributions.

        Parameters
        ----------
        squeeze : bool, optional
            Squeeze the summaries. If you only have one chain, squeeze will not return
            a length one list, just the single summary. If this is false, you will
            get a length one list.
        parameters : list[str], optional
            A list of parameters which to generate summaries for.
        chains : list[int|str], optional
            A list of the chains to get a summary of.

        Returns
        -------
        list of dictionaries
            One entry per chain, parameter bounds stored in dictionary with parameter as key
        """
        results = []
        if chains is None:
            chains = self.parent.get_mcmc_chains()
        else:
            if isinstance(chains, (int, str)):
                chains = [chains]
            if isinstance(chains[0], (int, str)):
                chains = [
                    self.parent.chains[i]
                    for c in chains
                    for i in self.parent._get_chain(c)
                ]

        for chain in chains:
            res = {}
            params_to_find = parameters if parameters is not None else chain.parameters
            for p in params_to_find:
                if p not in chain.parameters:
                    continue
                summary = self.get_parameter_summary(chain, p)
                res[p] = summary
            results.append(res)
        if squeeze and len(results) == 1:
            return results[0]
        return results

    def get_max_posteriors(self, parameters=None, squeeze=True, chains=None):
        """Gets the maximum posterior point in parameter space from the passed parameters.

        Requires the chains to have set `posterior` values.

        Parameters
        ----------
        parameters : str|list[str]
            The parameters to find
        squeeze : bool, optional
            Squeeze the summaries. If you only have one chain, squeeze will not return
            a length one list, just the single summary. If this is false, you will
            get a length one list.
        chains : list[int|str], optional
            A list of the chains to get a summary of.

        Returns
        -------
        list of two-tuples
            One entry per chain, two-tuple represents the max-likelihood coordinate
        """

        results = []
        if chains is None:
            chains = self.parent.chains
        else:
            if isinstance(chains, (int, str)):
                chains = [chains]
            chains = [
                self.parent.chains[i] for c in chains for i in self.parent._get_chain(c)
            ]

        if isinstance(parameters, str):
            parameters = [parameters]

        for chain in chains:
            if chain.posterior_max_index is None:
                results.append(None)
                continue
            res = {}
            params_to_find = parameters if parameters is not None else chain.parameters
            for p in params_to_find:
                if p in chain.parameters:
                    res[p] = chain.posterior_max_params[p]
            results.append(res)

        if squeeze and len(results) == 1:
            return results[0]
        return results

    def get_parameter_summary(self, chain, parameter):
        # Ensure config has been called so we get the statistics set in config
        if not self.parent._configured:
            self.parent.configure()
        callback = self._summaries[chain.config["statistics"]]
        return chain.get_summary(parameter, callback)

    def get_correlations(self, chain=0, parameters=None):
        """
        Takes a chain and returns the correlation between chain parameters.

        Parameters
        ----------
        chain : int|str, optional
            The chain index or name. Defaults to first chain.
        parameters : list[str], optional
            The list of parameters to compute correlations. Defaults to all parameters
            for the given chain.

        Returns
        -------
            tuple
                The first index giving a list of parameter names, the second index being the
                2D correlation matrix.
        """
        parameters, cov = self.get_covariance(chain=chain, parameters=parameters)
        diag = np.sqrt(np.diag(cov))
        divisor = diag[None, :] * diag[:, None]
        correlations = cov / divisor
        return parameters, correlations

    def get_covariance(self, chain=0, parameters=None):
        """
        Takes a chain and returns the covariance between chain parameters.

        Parameters
        ----------
        chain : int|str, optional
            The chain index or name. Defaults to first chain.
        parameters : list[str], optional
            The list of parameters to compute correlations. Defaults to all parameters
            for the given chain.

        Returns
        -------
            tuple
                The first index giving a list of parameter names, the second index being the
                2D covariance matrix.
        """
        index = self.parent._get_chain(chain)
        assert len(index) == 1, "Please specify only one chain, have %d chains" % len(
            index
        )
        chain = self.parent.chains[index[0]]
        if parameters is None:
            parameters = chain.parameters

        data = chain.get_data(parameters)
        cov = np.atleast_2d(np.cov(data, aweights=chain.weights, rowvar=False))

        return parameters, cov

    def get_correlation_table(
        self,
        chain=0,
        parameters=None,
        caption="Parameter Correlations",
        label="tab:parameter_correlations",
    ):
        """
        Gets a LaTeX table of parameter correlations.

        Parameters
        ----------
        chain : int|str, optional
            The chain index or name. Defaults to first chain.
        parameters : list[str], optional
            The list of parameters to compute correlations. Defaults to all parameters
            for the given chain.
        caption : str, optional
            The LaTeX table caption.
        label : str, optional
            The LaTeX table label.

        Returns
        -------
            str
                The LaTeX table ready to go!
        """
        parameters, cor = self.get_correlations(chain=chain, parameters=parameters)
        return self._get_2d_latex_table(parameters, cor, caption, label)

    def get_covariance_table(
        self,
        chain=0,
        parameters=None,
        caption="Parameter Covariance",
        label="tab:parameter_covariance",
    ):
        """
        Gets a LaTeX table of parameter covariance.

        Parameters
        ----------
        chain : int|str, optional
            The chain index or name. Defaults to first chain.
        parameters : list[str], optional
            The list of parameters to compute correlations. Defaults to all parameters
            for the given chain.
        caption : str, optional
            The LaTeX table caption.
        label : str, optional
            The LaTeX table label.

        Returns
        -------
            str
                The LaTeX table ready to go!
        """
        parameters, cov = self.get_covariance(chain=chain, parameters=parameters)
        return self._get_2d_latex_table(parameters, cov, caption, label)

    def _get_smoothed_histogram(self, chain, parameter):
        """Generate a smooth estimate of a 1D PDF from some samples using
        Kernel Density Estimation, correcting at the boundaries. (follows Joe
        Zuntz's cosmosis/density.py and getdist/mcsamples.py)
        """
        weights = chain.weights
        data = chain.get_data(parameter)
        neff = np.sum(weights) ** 2 / np.sum(weights**2)
        stdev = std_weight(data, weights)

        # This is a standard factor for 1D KDEs
        scott_factor = neff ** (-0.2)

        # We optionally allow user-defined additional smoothing to the optimal
        # smoothing
        width = 1.06 * stdev * scott_factor * neff ** (1.0 / 5 - 1.0 / (4 * 1 + 5))

        if chain.grid:
            edges = get_grid_bins(data)
            N = len(edges) - 1
        else:
            minv, maxv = get_extents(data, weights)
            N = 1024
            edges = np.linspace(minv, maxv, N + 1)

        x = 0.5 * (edges[1:] + edges[:-1])
        hist, _ = np.histogram(data, bins=edges, density=True, weights=weights)
        w = edges[1] - edges[0]
        # smoothing scale in units of the bin width
        s = width / w

        # Make the Gaussian kernel with which we are convolving.
        # We go out to 3 sigma
        window_width = min(max(1, int(2.5 * s)), N // 2 - 2)
        window_x = np.arange(-window_width, window_width + 1)
        kernel = np.exp(-0.5 * (window_x / s) ** 2)

        # Generate the smoothed version.  If we do not need
        # the boundary smoothing then this is our final output
        P_smooth = fftconvolve(hist, kernel, "same")
        P_final = P_smooth / P_smooth.sum() * (x[1] - x[0])

        correct_boundary = chain.config["correct_boundary"]
        if correct_boundary:
            # Generate the mask, a top-hat which cuts off where the
            # boundaries are
            full_width = N + 2 * window_width
            mask = np.ones(full_width)
            mask[:window_width] = 0
            mask[window_width] = 0.5
            mask[-window_width:] = 0
            mask[-(window_width + 1)] = 0.5

            # Linear boundary kernel correction Jones 1993 Jones and Foster 1996
            # implemented according to getdist/mcsamples.py
            a0 = fftconvolve(mask, kernel, "valid")

            # Avoid a divide-by-zero
            ix = np.nonzero(a0 * P_smooth)
            a0 = a0[ix]
            P_norm = P_smooth[ix] / a0
            xK = window_x * kernel
            x2K = xK * window_x
            a1 = fftconvolve(mask, xK, mode="valid")[ix]
            a2 = fftconvolve(mask, x2K, mode="valid")[ix]
            xP = fftconvolve(hist, xK, mode="same")[ix]

            # Apply the correction
            scaling = (P_smooth[ix] * a2 - xP * a1) / (a0 * a2 - a1**2)
            P_final = P_smooth.copy()
            P_final[ix] = P_norm * np.exp(np.minimum(scaling / P_norm, 4) - 1)

        correct_multbias = chain.config["correct_multbias"]
        if correct_multbias:
            mask = np.ones(N)
            mask[0] = 0.5
            mask[-1] = 0.5
            a0 = fftconvolve(mask, kernel, "same")

            tmp = P_final.copy()
            tmp[tmp == 0] = 1
            fine = hist / tmp
            conv = fftconvolve(fine, kernel, "same")
            P_final = P_final * conv / a0

        # Normalize and return
        P_final /= P_final.sum() * (x[1] - x[0])

        cs = P_final.cumsum()
        cs /= cs.max()
        return x, P_final, cs

    def _get_2d_latex_table(self, parameters, matrix, caption, label):
        latex_table = get_latex_table_frame(caption=caption, label=label)
        column_def = "c|%s" % ("c" * len(parameters))
        hline_text = "        \\hline\n"

        table = ""
        table += " & ".join([""] + parameters) + "\\\\ \n"
        table += hline_text
        max_len = max([len(s) for s in parameters])
        format_string = "        %%%ds" % max_len
        for p, row in zip(parameters, matrix):
            table += format_string % p
            for r in row:
                table += " & %5.2f" % r
            table += " \\\\ \n"
        table += hline_text
        return latex_table % (column_def, table)

    def get_parameter_text(self, lower, maximum, upper, wrap=False):
        """Generates LaTeX appropriate text from marginalised parameter bounds.

        Parameters
        ----------
        lower : float
            The lower bound on the parameter
        maximum : float
            The value of the parameter with maximum probability
        upper : float
            The upper bound on the parameter
        wrap : bool
            Wrap output text in dollar signs for LaTeX

        Returns
        -------
        str
            The formatted text given the parameter bounds
        """
        if lower is None or upper is None:
            return ""
        upper_error = upper - maximum
        lower_error = maximum - lower
        if upper_error != 0 and lower_error != 0:
            resolution = min(
                np.floor(np.log10(np.abs(upper_error))),
                np.floor(np.log10(np.abs(lower_error))),
            )
        elif upper_error == 0 and lower_error != 0:
            resolution = np.floor(np.log10(np.abs(lower_error)))
        elif upper_error != 0 and lower_error == 0:
            resolution = np.floor(np.log10(np.abs(upper_error)))
        else:
            resolution = np.floor(np.log10(np.abs(maximum)))
        factor = 0
        fmt = "%0.1f"
        r = 1
        if np.abs(resolution) > 2:
            factor = -resolution
        if resolution == 2:
            fmt = "%0.0f"
            factor = -1
            r = 0
        if resolution == 1:
            fmt = "%0.0f"
        if resolution == -1:
            fmt = "%0.2f"
            r = 2
        elif resolution == -2:
            fmt = "%0.3f"
            r = 3
        upper_error *= 10**factor
        lower_error *= 10**factor
        maximum *= 10**factor
        upper_error = round(upper_error, r)
        lower_error = round(lower_error, r)
        maximum = round(maximum, r)
        if maximum == -0.0:
            maximum = 0.0
        if resolution == 2:
            upper_error *= 10**-factor
            lower_error *= 10**-factor
            maximum *= 10**-factor
            factor = 0
            fmt = "%0.0f"
        upper_error_text = fmt % upper_error
        lower_error_text = fmt % lower_error
        if upper_error_text == lower_error_text:
            text = r"%s\pm %s" % (fmt, "%s") % (maximum, lower_error_text)
        else:
            text = (
                r"%s^{+%s}_{-%s}"
                % (fmt, "%s", "%s")
                % (maximum, upper_error_text, lower_error_text)
            )
        if factor != 0:
            text = r"\left( %s \right) \times 10^{%d}" % (text, -factor)
        if wrap:
            text = "$%s$" % text
        return text

    def get_parameter_summary_mean(self, chain, parameter):
        desired_area = chain.config["summary_area"]
        xs, _, cs = self._get_smoothed_histogram(chain, parameter)
        vals = [0.5 - desired_area / 2, 0.5, 0.5 + desired_area / 2]
        bounds = interp1d(cs, xs)(vals)
        bounds[1] = 0.5 * (bounds[0] + bounds[2])
        return bounds

    def get_parameter_summary_cumulative(self, chain, parameter):
        xs, _, cs = self._get_smoothed_histogram(chain, parameter)
        desired_area = chain.config["summary_area"]
        vals = [0.5 - desired_area / 2, 0.5, 0.5 + desired_area / 2]
        bounds = interp1d(cs, xs)(vals)
        return bounds

    def get_parameter_summary_max(self, chain, parameter):
        xs, ys, cs = self._get_smoothed_histogram(chain, parameter)
        desired_area = chain.config["summary_area"]
        n_pad = 1000
        x_start = xs[0] * np.ones(n_pad)
        x_end = xs[-1] * np.ones(n_pad)
        y_start = np.linspace(0, ys[0], n_pad)
        y_end = np.linspace(ys[-1], 0, n_pad)
        xs = np.concatenate((x_start, xs, x_end))
        ys = np.concatenate((y_start, ys, y_end))
        cs = ys.cumsum()
        cs = cs / cs.max()
        startIndex = ys.argmax()
        maxVal = ys[startIndex]
        minVal = 0
        threshold = 0.003
        x1 = None
        x2 = None
        count = 0
        while x1 is None:
            mid = (maxVal + minVal) / 2.0
            count += 1
            try:
                if count > 50:
                    raise ValueError("Failed to converge")
                i1 = startIndex - np.where(ys[:startIndex][::-1] < mid)[0][0]
                i2 = startIndex + np.where(ys[startIndex:] < mid)[0][0]
                area = cs[i2] - cs[i1]
                deviation = np.abs(area - desired_area)
                if deviation < threshold:
                    x1 = xs[i1]
                    x2 = xs[i2]
                elif area < desired_area:
                    maxVal = mid
                elif area > desired_area:
                    minVal = mid
            except ValueError:
                self._logger.warning(
                    "Parameter %s in chain %s is not constrained"
                    % (parameter, chain.name)
                )
                return [None, xs[startIndex], None]

        return [x1, xs[startIndex], x2]

    def get_paramater_summary_max_symmetric(self, chain, parameter):
        xs, ys, cs = self._get_smoothed_histogram(chain, parameter)
        desired_area = chain.config["summary_area"]

        x_to_c = interp1d(xs, cs, bounds_error=False, fill_value=(0, 1))

        # Get max likelihood x
        max_index = ys.argmax()
        x = xs[max_index]

        # Estimate width
        h = 0.5 * (xs[-1] - xs[0])
        prev_h = 0

        # Hone in on right answer
        while True:
            current_area = x_to_c(x + h) - x_to_c(x - h)
            if np.abs(current_area - desired_area) < 0.0001:
                return [x - h, x, x + h]
            temp = h
            h += 0.5 * np.abs(prev_h - h) * (1 if current_area < desired_area else -1)
            prev_h = temp

    def get_parameter_summary_max_shortest(self, chain, parameter):
        xs, ys, cs = self._get_smoothed_histogram(chain, parameter)
        desired_area = chain.config["summary_area"]

        c_to_x = interp1d(cs, xs, bounds_error=False, fill_value=(-np.inf, np.inf))

        # Get max likelihood x
        max_index = ys.argmax()
        x = xs[max_index]

        # Pair each lower bound with an upper to get the right area
        x2 = c_to_x(cs + desired_area)
        dists = x2 - xs
        mask = (xs > x) | (x2 < x)  # Ensure max point is inside the area
        dists[mask] = np.inf
        ind = dists.argmin()
        return [xs[ind], x, x2[ind]]

    def get_parameter_summary_max_central(self, chain, parameter):
        xs, ys, cs = self._get_smoothed_histogram(chain, parameter)
        desired_area = chain.config["summary_area"]

        c_to_x = interp1d(cs, xs)

        # Get max likelihood x
        max_index = ys.argmax()
        x = xs[max_index]

        vals = [0.5 - 0.5 * desired_area, 0.5 + 0.5 * desired_area]
        xvals = c_to_x(vals)

        return [xvals[0], x, xvals[1]]
