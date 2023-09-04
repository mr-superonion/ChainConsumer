# -*- coding: utf-8 -*-
from matplotlib.colors import rgb2hex
import matplotlib.pyplot as plt
import numpy as np

# Colours drawn from material designs colour pallet at https://material.io/guidelines/style/color.html


class Colors(object):
    def __init__(self):
        self.color_map = {
            "blue": "#252c8d",  # blue
            "red": "#ac171f",  # red
            "green": "#0e7003",  # green
            "orange": "#ed611e",  # orange
            "purple": "#9700ff",  # purple
            "brown": "#ae7b26",  # brown
            "pink": "#ff7d9c",  # pink
            "magenta": "#c01e62",  # magenta
            "grey": "#3d3d3d",  # grey
            "yellow": "#fecc2e",  # yellow
            "teal": "#52bcbc",  # teal
            "oliver": "#7f7f00",  # olive
        }
        self.aliases = {
            "k": "black",
            "b": "blue",
            "r": "red",
            "g": "green",
            "o": "orange",
            "p": "purple",
            "m": "magenta",
            "e": "grey",
            "y": "yellow",
            "t": "teal",
            "br": "brown",
            "pk": "pink",
            "ol": "oliver",
        }
        self.default_colors = [
            "blue",
            "red",
            "green",
            "orange",
            "purple",
            "brown",
            "pink",
            "magenta",
            "grey",
            "yellow",
            "teal",
            "oliver",
        ]

    def format(self, color):
        if isinstance(color, np.ndarray):
            color = rgb2hex(color)
        if color[0] == "#":
            return color
        elif color in self.color_map:
            return self.color_map[color]
        elif color in self.aliases:
            alias = self.aliases[color]
            return self.color_map[alias]
        else:
            raise ValueError("Color %s is not mapped. Please give a hex code" % color)

    def get_formatted(self, list_colors):
        return [self.format(c) for c in list_colors]

    def get_default(self):
        return self.get_formatted(self.default_colors)

    def get_colormap(self, num, cmap_name, scale=0.7):  # pragma: no cover
        color_list = self.get_formatted(
            plt.get_cmap(cmap_name)(np.linspace(0.05, 0.9, num))
        )
        scales = scale + (1 - scale) * np.abs(1 - np.linspace(0, 2, num))
        scaled = [self.scale_colour(c, s) for c, s in zip(color_list, scales)]
        return scaled

    def scale_colour(self, colour, scalefactor):  # pragma: no cover
        if isinstance(colour, np.ndarray):
            r, g, b = colour[:3] * 255.0
        else:
            hexx = colour.strip("#")
            if scalefactor < 0 or len(hexx) != 6:
                return hexx
            r, g, b = int(hexx[:2], 16), int(hexx[2:4], 16), int(hexx[4:], 16)
        r = self._clamp(int(r * scalefactor))
        g = self._clamp(int(g * scalefactor))
        b = self._clamp(int(b * scalefactor))
        return "#%02x%02x%02x" % (r, g, b)

    def _clamp(self, val, minimum=0, maximum=255):
        if val < minimum:
            return minimum
        if val > maximum:
            return maximum
        return val
