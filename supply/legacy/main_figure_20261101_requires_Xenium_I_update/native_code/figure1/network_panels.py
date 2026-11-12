"""Vector network panels with shared node and edge scales for both methods."""
import numpy as np
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path


def draw_network(ax, groups, weights, max_weight):
    names = list(groups.group)
    theta = np.arange(len(names)) * 2 * np.pi / len(names)
    positions = {name: np.array([np.cos(t), np.sin(t)]) for name, t in zip(names, theta)}
    colors = dict(zip(groups.group, groups.colour))
    counts = dict(zip(groups.group, groups.cells))
    for row in weights.sort_values('weight').itertuples():
        start, end = positions[row.sender], positions[row.receiver]
        width = .16 + 1.25 * row.weight / max_weight
        options = dict(arrowstyle='-|>', mutation_scale=3.4,
                       linewidth=width, color=colors[row.sender], alpha=.62, zorder=1)
        if row.sender == row.receiver:
            tangent = np.array([-start[1], start[0]])
            path = Path([start + .05 * tangent,
                         start + .38 * start + .27 * tangent,
                         start + .38 * start - .27 * tangent,
                         start - .05 * tangent],
                        [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
            arrow = FancyArrowPatch(path=path, **options)
        else:
            arrow = FancyArrowPatch(start, end, connectionstyle='arc3,rad=.10',
                                    shrinkA=2.5, shrinkB=2.5, **options)
        ax.add_patch(arrow)
    for name in names:
        x, y = positions[name]
        ax.scatter([x], [y], s=8 + 37 * counts[name] / max(counts.values()),
                   c=[colors[name]], linewidths=.25, edgecolors='white', zorder=3)
        label = 1.13 * positions[name]
        ax.text(*label, name, fontsize=6.0,
                ha='left' if x > .15 else ('right' if x < -.15 else 'center'),
                va='bottom' if y > .8 else ('top' if y < -.8 else 'center'),
                zorder=4)
    ax.set(xlim=(-1.88, 1.88), ylim=(-1.57, 1.57), aspect='equal')
    ax.axis('off')
