"""New figures: Fig. 1 (capture problem and model), Fig. 3 (computed shapes and energy split) and the graphical
abstract. Inputs: datos/fem_shapes_channel.json, datos/fem_shapes_channel_profile.json,
../../webapp/energy_curves.json. Output: ../../figures/ at 600 dpi."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
from matplotlib import cm, colors

AQUI = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(AQUI, "..", "datos")
OUT = os.path.abspath(os.path.join(AQUI, "..", "..", "figures")); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7, "axes.linewidth": 0.6, "savefig.dpi": 600,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6})
MM = 1 / 25.4
C_CORTEX, C_CYT, C_ENV, C_NUC = "#c2416b", "#f2c6d3", "#3b2f70", "#7b6cc4"
C_LEU, C_TUM, C_WALL = "#7b6cc4", "#c2416b", "#a7aeb8"
SH = json.load(open(os.path.join(D, "fem_shapes_channel.json")))
PROF = json.load(open(os.path.join(D, "fem_shapes_channel_profile.json")))["profile"]
SCALE = lambda s: np.minimum(s / (s + 1), 5 / 6) * 1.2          # strain 0, .5, 1, 2, >=5 -> 0, .4, .6, .8, 1


def Rw(z, Rn, p=PROF):
    g = lambda z: 1 + 0.5 * (np.tanh((z - p["zt"]) / p["L"]) - np.tanh((z + p["zt"]) / p["L"]))
    return Rn + (p["R_ext"] - Rn) * (g(np.abs(z)) - g(0)) / (1 - g(0))


def letter(ax, s, x=-0.02, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def cell_polys(ref, frame, zc=0.0):
    P = np.array(ref["mesh"]["p"]); T = np.array(ref["mesh"]["t"]); a = ref["a"]
    u = np.zeros_like(P) if frame is None else np.array(frame["u"])
    xr, xz = (P[0] + u[0]) * a, (P[1] + u[1]) * a
    polys, idx = [], []
    for sr in (1, -1):
        for sz in (1, -1):
            for e, t in enumerate(T.T):
                polys.append(np.column_stack([zc + sz * xz[t], sr * xr[t]])); idx.append(e)
    return polys, np.array(idx)


def wall(ax, Rn, span, color=C_WALL):
    z = np.linspace(-span, span, 800)
    for sg in (1, -1):
        ax.fill_between(z, sg * Rw(z, Rn), sg * (PROF["R_ext"] + 3), color=color, lw=0)
        ax.plot(z, sg * Rw(z, Rn), color="#5d6672", lw=0.6)


# =============================================================== Fig. 1: capture problem and model
fig = plt.figure(figsize=(174 * MM, 120 * MM))
gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.25, 1.1], height_ratios=[0.85, 1], wspace=0.42, hspace=0.32)

# (a) the capture problem
ax = fig.add_subplot(gs[0, :]); letter(ax, "a", y=1.0); ax.set_xlim(0, 100); ax.set_ylim(-21, 14); ax.axis("off")
z = np.linspace(0, 100, 900)
prof_a = lambda x: 9.5 - 4.6 * np.exp(-((x - 62) / 6.5) ** 2)
for sg in (1, -1):
    ax.fill_between(z, sg * prof_a(z), sg * 14, color=C_WALL, lw=0); ax.plot(z, sg * prof_a(z), color="#5d6672", lw=0.6)
rng = np.random.default_rng(3)
for x, y in [(6, 3.5), (13, -4.2), (38, -3.8), (44, 4.2), (51, -1.5), (77, 0.8), (85, -3.5), (93, 3)]:
    ax.add_patch(Circle((x, y), 2.0, fc="#ddd8f3", ec=C_LEU, lw=0.6)); ax.add_patch(Circle((x + 0.3, y), 1.15, fc=C_LEU, ec="none"))
ax.add_patch(Circle((27, 0.4), 4.6, fc=C_CYT, ec=C_CORTEX, lw=1.4)); ax.add_patch(Circle((27.3, 0.4), 2.3, fc=C_NUC, ec=C_ENV, lw=0.7))
ax.text(27, -6.4, "CTC", ha="center", fontsize=6.5, color=C_TUM)
ax.text(8, -7.6, "leucocytes", ha="center", fontsize=6.5, color=C_LEU)
ax.add_patch(FancyArrowPatch((2, 11.5), (16, 11.5), arrowstyle="-|>", mutation_scale=8, color="#4d5868", lw=0.8))
ax.text(9, 12.4, "flow", ha="center", fontsize=6.5, color="#4d5868")
ax.text(62, 11.3, "constriction or pore", ha="center", fontsize=6.5, color="#4d5868")
qs = [(17, "1  Size\nwhich cells touch the wall?"), (50, "2  Cell\nwhich compartment resists entry?"),
      (83, "3  Device\nhow do geometry and flow set the pressure?")]
for x, s in qs:
    ax.text(x, -17.5, s, ha="center", va="center", fontsize=6.6, color="#1d2430", linespacing=1.3,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#c9ced6", lw=0.6))

# (b) four-compartment cell
ax = fig.add_subplot(gs[1, 0]); letter(ax, "b"); ax.set_aspect("equal"); ax.axis("off")
ref = SH["tumour"]; polys, idx = cell_polys(ref, None)
nuc = np.array(ref["mesh"]["nuc"], bool)
ax.add_collection(PolyCollection(polys, facecolors=[C_NUC if nuc[e] else C_CYT for e in idx], edgecolors="face", lw=0.2))
th = np.linspace(0, 2 * np.pi, 300); a = ref["a"]; an = ref["an_a"] * a
ax.plot(a * np.cos(th), a * np.sin(th), color=C_CORTEX, lw=2.2); ax.plot(an * np.cos(th), an * np.sin(th), color=C_ENV, lw=1.2)
ax.plot([0, 0], [-a * 1.15, a * 1.15], color="#4d5868", lw=0.5, ls=(0, (3, 2)))
lab = [((a * 0.70, a * 0.72), (a * 1.18, a * 1.12), "cortex\n(Skalak shell)", C_CORTEX),
       ((a * 0.62, -a * 0.38), (a * 1.18, -a * 0.75), "cytoplasm\n(neo-Hookean)", "#8c2d4f"),
       ((-an * 0.72, an * 0.70), (-a * 1.25, a * 1.22), "nuclear\nenvelope", C_ENV),
       ((-an * 0.3, -an * 0.4), (-a * 1.25, -a * 1.32), "nucleoplasm", C_ENV)]
for xy, xyt, s, c in lab:
    ax.annotate(s, xy=xy, xytext=xyt, fontsize=6, color=c, va="center", linespacing=1.15,
                arrowprops=dict(arrowstyle="-", color=c, lw=0.5))
ax.text(-a * 0.06, a * 1.32, "symmetry axis", fontsize=5.3, color="#4d5868", ha="right", va="center")
ax.set_xlim(-a * 1.3, a * 2.15); ax.set_ylim(-a * 1.55, a * 1.45)

# (c) device geometry
ax = fig.add_subplot(gs[1, 1]); letter(ax, "c")
a = 8.74; Rt = a / 1.70; Rwide = 1.6 * a; zt, L = 4.0, 5.0
g0 = 0.5 * (1 + np.tanh(-zt / L)); Rmin = Rt + g0 * (Rwide - Rt)
zz = np.linspace(-22, 22, 600)
Rwz = Rt + (Rwide - Rt) * 0.5 * (1 + np.tanh((np.abs(zz) - zt) / L))
ax.fill_between(zz, Rwz, Rwide + 3, color=C_WALL, lw=0); ax.plot(zz, Rwz, color="#5d6672", lw=0.7)
ax.plot(a * np.cos(th), a * np.sin(th) * 0 + a * np.abs(np.sin(th)), color=C_TUM, lw=0.8, ls=(0, (3, 2)))
ax.axhline(0, color="#4d5868", lw=0.5, ls=(0, (3, 2)))
for y, dy, s, c in [(Rt, -0.9, "$R_\\mathrm{t}$ nominal, not reached", "#4d5868"), (Rmin, 0.5, "$R_\\mathrm{min}$ at the neck", "#1d2430"), (Rwide, 0.5, "$R_\\mathrm{wide}=1.6a$", "#4d5868")]:
    ax.plot([-22, 22], [y, y], color=c, lw=0.4, ls=":"); ax.text(-21.5, y + dy, s, fontsize=5.6, va="center", color=c)
ax.annotate("", xy=(-zt, 1.2), xytext=(zt, 1.2), arrowprops=dict(arrowstyle="<->", lw=0.6)); ax.text(0, 1.9, "$2z_\\mathrm{t}$", ha="center", fontsize=6)
ax.annotate("", xy=(zt, 3.4), xytext=(zt + L, 3.4), arrowprops=dict(arrowstyle="<->", lw=0.6)); ax.text(zt + L / 2, 4.1, "$L_\\mathrm{tap}$", ha="center", fontsize=6)
ax.text(-a * 0.8, a * 0.62, "cell, radius $a$", fontsize=5.6, color=C_TUM, ha="right")
ax.set_xlim(-22, 22); ax.set_ylim(0, Rwide + 3); ax.set_xlabel("axial position $z$ (µm)"); ax.set_ylabel("radius (µm)")
ax.spines[["top", "right"]].set_visible(False)

# (d) from energy to entry pressure
ax = fig.add_subplot(gs[1, 2]); letter(ax, "d")
EC = json.load(open(os.path.join(AQUI, "..", "..", "webapp", "energy_curves.json")))["channel_PDMS"]
from scipy.interpolate import PchipInterpolator
u = np.array(EC["u"]); gE = np.array(EC["g"]); o = np.argsort(u); u, gE = u[o], gE[o]
keep = u < 1; u, gE = np.append(u[keep], 1.0), np.append(gE[keep], 0.0); gE = np.maximum.accumulate(gE[::-1])[::-1]
sp = PchipInterpolator(u, gE); pr = EC["profile"]; Rn = a / 1.45
zz = np.linspace(-14, 0.0, 1400); Rz = Rw(zz, Rn, pr); uz = Rz / a; contact = uz < 1
W = np.where(contact, sp(np.clip(uz, u.min(), 1)), 0.0)
dRdz = np.gradient(Rz, zz); q = np.where(contact, np.abs(sp.derivative()(np.clip(uz, u.min(), 1))) * np.abs(dRdz) / (np.pi * uz ** 2), 0.0)
ax.plot(zz, W / W.max(), color="#1d2430", lw=1.1, label="stored energy $W(s)$")
ax.plot(zz, q / q.max(), color=C_TUM, lw=1.1, label="resisting pressure")
i = int(np.argmax(q)); ax.plot(zz[i], 1, "o", color=C_TUM, ms=3.5)
ax.annotate("$\\Delta P_\\mathrm{crit}$ = maximum\nalong the path", xy=(zz[i], 1), xytext=(-13.5, 0.66), fontsize=5.8, color=C_TUM,
            arrowprops=dict(arrowstyle="-", color=C_TUM, lw=0.5))
ax.set_xlabel("position of the cell centre $s$ (µm)"); ax.set_ylabel("normalised value"); ax.set_ylim(0, 1.12)
ax.legend(fontsize=5.6, frameon=False, loc="upper left"); ax.spines[["top", "right"]].set_visible(False)
ax.text(-0.3, 0.04, "neck", ha="right", fontsize=5.6, color="#4d5868")
fig.savefig(os.path.join(OUT, "Fig1.png"), bbox_inches="tight"); plt.close(fig)

# =============================================================== Fig. 3: computed shapes and energy split
fig = plt.figure(figsize=(174 * MM, 92 * MM))
gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 0.95], wspace=0.08, hspace=0.42)
norm = colors.Normalize(0, 1); cmap = matplotlib.colormaps["cividis"]
for row, (name, lab) in enumerate([("tumour", "breast tumour cell (a = 8.74 µm)"), ("leucocyte", "leucocyte (a = 4.20 µm)")]):
    ref = SH[name]; a = ref["a"]
    for col, fr in enumerate(ref["frames"]):
        ax = fig.add_subplot(gs[row, col]); ax.set_aspect("equal"); ax.axis("off")
        Rn = a / fr["conf"]; span = 2.3 * a
        wall(ax, Rn, span)
        polys, idx = cell_polys(ref, fr)
        s = np.array(fr["strain"])[idx]
        ax.add_collection(PolyCollection(polys, array=SCALE(s), cmap=cmap, norm=norm, edgecolors="face", lw=0.15))
        ax.set_xlim(-span, span); ax.set_ylim(-1.75 * a, 1.75 * a)
        if row == 0:
            ax.set_title(f"$a/R_\\mathrm{{min}}$ = {fr['conf']:.2f}", fontsize=6.8, pad=3)
        if col == 0:
            ax.text(-span * 1.04, 0, lab.replace(" (", "\n("), fontsize=6.3, fontweight="bold", va="center", ha="right", rotation=90)
            letter(ax, "ab"[row], x=-0.2, y=1.0)
    ax = fig.add_subplot(gs[row, 3]); letter(ax, "cd"[row], x=-0.18, y=1.05)
    comp = [("W_cortex", "cortex", C_CORTEX), ("W_cytoplasm", "cytoplasm", C_CYT), ("W_envelope", "nuclear envelope", C_ENV), ("W_nucleus", "nucleoplasm", C_NUC)]
    left = np.zeros(len(ref["frames"]))
    for k, nm, c in comp:
        v = np.array([f["Wc"][k] / f["W"] for f in ref["frames"]])
        ax.barh(range(len(v)), v, left=left, color=c, height=0.62, label=nm, edgecolor="white", lw=0.4)
        left += v
    ax.set_yticks(range(len(left))); ax.set_yticklabels([f"{f['conf']:.2f}" for f in ref["frames"]], fontsize=6)
    ax.set_xlim(0, 1); ax.invert_yaxis(); ax.set_xlabel("share of stored energy", fontsize=6.5)
    ax.set_ylabel("$a/R_\\mathrm{min}$", fontsize=6.5, labelpad=2); ax.spines[["top", "right"]].set_visible(False)
    for j, f in enumerate(ref["frames"]):
        ax.text(0.02, j, f"{f['Wc']['W_cortex'] / f['W']:.2f}", va="center", fontsize=5.6, color="white")
    if row == 0:
        ax.legend(fontsize=5.4, frameon=False, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.08), handlelength=1)
cax = fig.add_axes([0.13, 0.04, 0.42, 0.018])
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cb.set_ticks([SCALE(np.array(v)) for v in (0, 0.5, 1, 2, 5)]); cb.set_ticklabels(["0", "0.5", "1", "2", "≥5"])
cb.set_label("equivalent Green–Lagrange strain", fontsize=6.5); cb.ax.tick_params(labelsize=6)
fig.savefig(os.path.join(OUT, "Fig3.png"), bbox_inches="tight"); plt.close(fig)

# =============================================================== graphical abstract (1328 x 531 px)
fig = plt.figure(figsize=(1328 / 300, 531 / 300), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis("off")
ref = SH["tumour"]; fr = ref["frames"][-1]; a = ref["a"]; Rn = a / fr["conf"]
sub = fig.add_axes([0.01, 0.12, 0.36, 0.78]); sub.set_aspect("equal"); sub.axis("off")
wall(sub, Rn, 2.2 * a); polys, idx = cell_polys(ref, fr)
sub.add_collection(PolyCollection(polys, array=SCALE(np.array(fr["strain"])[idx]), cmap=cmap, norm=norm, edgecolors="face", lw=0.1))
sub.set_xlim(-2.2 * a, 2.2 * a); sub.set_ylim(-1.7 * a, 1.7 * a)
ax.text(19, 37, "computed cell at the neck", ha="center", fontsize=6.5, color="#4d5868")
msgs = [(39, 27, "Size is not enough", "patient-sized breast CTCs:\nonly 1.3–1.7-fold enrichment\nby size"),
        (39, 9, "The cortex resists entry", "while the aperture is wider\nthan the nucleus"),
        (69.5, 27, "Retain the tumour cell", "in these data tumour cells\nare stiffer than leucocytes:\nΔP$_\\mathrm{leu}$ < ΔP$_\\mathrm{applied}$ < ΔP$_\\mathrm{tum}$"),
        (69.5, 9, "Design by geometry and β", "entry shape and pressure\nper unit velocity set\nthe operating point")]
for x, y, h, txt in msgs:
    ax.text(x, y + 3.0, h, fontsize=6.3, fontweight="bold", color="#1d2430", va="bottom")
    ax.text(x, y + 2.4, txt, fontsize=5.2, color="#4d5868", va="top", linespacing=1.25)
fig.savefig(os.path.join(OUT, "graphical_abstract.png"), dpi=300); plt.close(fig)
print("Fig1.png, Fig3.png, graphical_abstract.png ->", OUT)
