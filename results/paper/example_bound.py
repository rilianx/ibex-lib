"""Illustration of the Lagrangian bound on a 1-D example (figure for the paper).

    min f(x) = x^4 - 3x^2 + x   s.t.  g(x) = -x - 1 <= 0,   x in [-2, 2]

The constrained optimum is x* = -1 (g active), f* = -3; KKT: f'(-1) + mu g'(-1) = 0
gives mu = 3, so L(x) = f(x) + 3 g(x) = x^4 - 3x^2 - 2x - 3, L'' = 12x^2 - 6.
The local solver returns (x^, mu) = (-1, 3); the incumbent is loup = f* = -3.
On a node box [a, b], with x_bar the projection of x^ onto it and
d = max(0, -min_{[a,b]} L'') (1-D Gerschgorin), the bound is
    beta = L(x_bar) + min_{[a,b]} L'(x_bar)(x - x_bar) - d/2 max_{[a,b]} (x - x_bar)^2.
Plain floating point here (the solver uses interval arithmetic).
Usage: python3 results/paper/example_bound.py [out.png]"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else 'results/paper/example-bound.png'
f = lambda x: x**4 - 3*x**2 + x
mu, xh, loup = 3.0, -1.0, -3.0
L = lambda x: f(x) + mu*(-x - 1)
dL = lambda x: 4*x**3 - 6*x + 1 - mu
d2L = lambda x: 12*x**2 - 6

def bound(a, b):
    xb = min(max(xh, a), b)
    h_lo = d2L(min(max(0.0, a), b))          # L'' = 12x^2 - 6 is smallest at the point closest to 0
    d = max(0.0, -h_lo)
    g = dL(xb)
    lin = min(g*(a - xb), g*(b - xb))
    pen = 0.5*d*max((a - xb)**2, (b - xb)**2)
    return xb, d, L(xb) + lin - pen

BLUE, ORANGE, INK, MUTED, GRID, SHADE = '#2a78d6', '#eb6834', '#1f1f1e', '#6b6a64', '#e6e5e0', '#f1f0ec'
X = np.linspace(-2, 2, 801)

def base(ax, title, note, xlim=(-2, 2), ylim=(-8.5, 3), loup_x=None, above=False):
    ax.axvspan(xlim[0], -1, color=SHADE, lw=0)
    ax.text(-1 - 0.012*(xlim[1]-xlim[0]), ylim[1] - 0.03*(ylim[1]-ylim[0]), 'infeasible (g > 0)',
            fontsize=7, color=MUTED, rotation=90, ha='right', va='top')
    ax.plot(X, f(X), color=INK, lw=1.6)
    ax.plot(X, L(X), color=BLUE, lw=2)
    ax.axhline(loup, color=MUTED, lw=1, ls=':')
    lx = loup_x if loup_x is not None else xlim[1]
    dy = 0.012*(ylim[1]-ylim[0])
    ax.text(lx, loup + dy if above else loup - dy, 'loup = f* = \u22123', fontsize=7.5, color=MUTED,
            ha='right', va='bottom' if above else 'top')
    ax.plot([xh], [L(xh)], 'o', color=INK, ms=6, zorder=5)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_title(title + '\n' + note, fontsize=8.5, color=INK, loc='left', linespacing=1.5)
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s_ in ('top', 'right'): ax.spines[s_].set_visible(False)
    for s_ in ('left', 'bottom'): ax.spines[s_].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=7.5)

def node(ax, a, b):
    xb, d, beta = bound(a, b)
    ax.axvspan(a, b, color=ORANGE, alpha=0.10, lw=0)
    xs = np.linspace(a, b, 200)
    ax.plot(xs, L(xb) + dL(xb)*(xs - xb) - 0.5*d*(xs - xb)**2, color=ORANGE, lw=2)
    ax.plot([a, b], [beta, beta], color=ORANGE, lw=1.4, ls='--')
    ax.annotate('\u03b2 = %.2f' % beta, (a, beta), xytext=(2, 3), textcoords='offset points',
                fontsize=8, color=INK, va='bottom')
    if xb != xh:
        ax.plot([xb], [L(xb)], 's', color=ORANGE, ms=5, zorder=5)
    return d, beta

fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.0))
ax = axes[0][0]
base(ax, '(a) Problem; local solution x\u0302 = \u22121 with multiplier \u03bc = 3',
     "L \u2264 f on the feasible set (weak duality); L(x\u0302) = f*, L'(x\u0302) = 0 (KKT)", loup_x=1.0, above=True)
ax.annotate('x\u0302', (xh, L(xh)), xytext=(4, 6), textcoords='offset points', fontsize=9, color=INK)

ax = axes[0][1]
base(ax, '(b) Box [\u22121.3, \u22120.8] around x\u0302: L convex, certified',
     "L'' \u2265 1.7 on the box, so d = 0: the bound is exact, \u03b2 = f*",
     xlim=(-1.6, -0.5), ylim=(-3.6, -1.5), loup_x=-0.52)
node(ax, -1.3, -0.8)

ax = axes[1][0]
d, beta = bound(-1.2, 0.2)[1:]
base(ax, '(c) Wide box [\u22121.2, 0.2]: L not convex, diagonal shift',
     "min L'' = \u22126 on the box, so d = %.0f: valid but weak, \u03b2 = %.2f < f*" % (d, beta), loup_x=1.95)
node(ax, -1.2, 0.2)

ax = axes[1][1]
d, beta = bound(-0.6, -0.4)[1:]
base(ax, '(d) Box [\u22120.6, \u22120.4] not containing x\u0302: pruned',
     "x\u0304 = projection of x\u0302 (square), d = %.2f: \u03b2 = %.2f > loup" % (d, beta),
     xlim=(-1.1, -0.2), ylim=(-3.6, -2.2), loup_x=-0.22)
node(ax, -0.6, -0.4)

for a_ in axes[1]: a_.set_xlabel('x', fontsize=8.5, color=MUTED)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
H = [Line2D([], [], color=INK, lw=1.6, label='f'),
     Line2D([], [], color=BLUE, lw=2, label='Lagrangian L = f + \u03bc g'),
     Line2D([], [], color=ORANGE, lw=2, label="tangent of L at x\u0304 \u2212 (d/2)(x \u2212 x\u0304)\u00b2"),
     Line2D([], [], color=ORANGE, lw=1.4, ls='--', label='bound \u03b2 on the box'),
     Patch(color=ORANGE, alpha=0.15, lw=0, label='node box')]
fig.legend(handles=H, loc='upper center', ncol=5, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.955))
fig.savefig(OUT, dpi=200)
fig.savefig(OUT.replace('.png', '.pdf'))
for a, b in ((-1.3, -0.8), (-1.2, 0.2), (-0.6, -0.4)):
    print('box [%g, %g]: x_bar %g d %.3f beta %.4f' % ((a, b) + bound(a, b)))
print('wrote', OUT)
