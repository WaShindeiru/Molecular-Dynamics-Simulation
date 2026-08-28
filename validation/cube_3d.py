import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

NX, NY, NZ = 8, 5, 4
X_SPLITS, Y_SPLITS = 3, 2

color_id = np.full((NX, NY, NZ), -1, dtype=int)

# x: [0,3) [3,6) [6,8)   y: [0,3) [3,5)   z: [0,2) [2,4)
color_id[0:3, 0:3, 0:2] = 0
color_id[3:6, 0:3, 0:2] = 1
color_id[6:8, 0:3, 0:2] = 2
color_id[0:3, 3:5, 0:2] = 3
color_id[3:6, 3:5, 0:2] = 4
color_id[6:8, 3:5, 0:2] = 5
color_id[0:3, 0:3, 2:4] = 6
color_id[3:6, 0:3, 2:4] = 7
color_id[6:8, 0:3, 2:4] = 8
color_id[0:3, 3:5, 2:4] = 9
color_id[3:6, 3:5, 2:4] = 10
color_id[6:8, 3:5, 2:4] = 11

PALETTE = [
    "#4C78A8", "#F58518", "#E45756", "#72B7B2",
    "#54A24B", "#EECA3B", "#B279A2", "#FF9DA6",
    "#9D755D", "#BAB0AC", "#1F77B4", "#D62728",
]

# placeholder sizes / spacing — tweak these
LABEL_FONTSIZE = 18
TICK_FONTSIZE = 15
TITLE_FONTSIZE = 20
LABEL_PAD = 10
TITLE_PAD_2D = 16
# 3D title: pad is ignored by mplot3d. Use axes fraction instead (1.0 = top).
TITLE_Y_3D = 0.97
# 3D "z" in axes fractions (0–1). x: toward the numbered z-axis, y: height along it.
Z_LABEL_AXES_X = 1.1
Z_LABEL_AXES_Y = 0.45


def apply_font(ax, title_pad=None, title_y=None):
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.xaxis.label.set_size(LABEL_FONTSIZE)
    ax.yaxis.label.set_size(LABEL_FONTSIZE)
    ax.xaxis.labelpad = LABEL_PAD
    ax.yaxis.labelpad = LABEL_PAD
    if hasattr(ax, "zaxis"):
        ax.zaxis.label.set_size(LABEL_FONTSIZE)
        ax.zaxis.labelpad = LABEL_PAD
        ax.tick_params(axis="z", labelsize=TICK_FONTSIZE)
    title = ax.get_title()
    if title:
        kwargs = {"fontsize": TITLE_FONTSIZE}
        if title_y is not None:
            kwargs["y"] = title_y
        elif title_pad is not None:
            kwargs["pad"] = title_pad
        ax.set_title(title, **kwargs)


def cube_faces(x, y, z, s):
    p = np.array([
        [x, y, z], [x + s, y, z], [x + s, y + s, z], [x, y + s, z],
        [x, y, z + s], [x + s, y, z + s], [x + s, y + s, z + s], [x, y + s, z + s],
    ])
    return [
        p[[0, 1, 2, 3]], p[[4, 5, 6, 7]],
        p[[0, 1, 5, 4]], p[[2, 3, 7, 6]],
        p[[1, 2, 6, 5]], p[[0, 3, 7, 4]],
    ]


def plot_3d(ax, color_id, gap=0.18, alpha=0.92):
    nx, ny, nz = color_id.shape
    s, off = 1.0 - gap, gap / 2.0
    used = sorted(int(v) for v in np.unique(color_id) if v >= 0)
    for cid in used:
        faces = []
        xs, ys, zs = np.where(color_id == cid)
        for x, y, z in zip(xs, ys, zs):
            faces.extend(cube_faces(x + off, y + off, z + off, s))
        ax.add_collection3d(Poly3DCollection(
            faces,
            facecolors=to_rgba(PALETTE[cid % len(PALETTE)], alpha=alpha),
            edgecolors="k",
            linewidths=0.35,
        ))
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    ax.set_zlim(0, nz)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("")
    ax.set_title("Podział zadań - widok 3d")
    ax.set_box_aspect((nx, ny, nz))
    ax.view_init(elev=26, azim=-40)
    apply_font(ax, title_y=TITLE_Y_3D)
    ax.text2D(
        Z_LABEL_AXES_X,
        Z_LABEL_AXES_Y,
        "z",
        transform=ax.transAxes,
        fontsize=LABEL_FONTSIZE,
        ha="center",
        va="center",
    )


def plot_top(ax, color_id, gap=0.18, z_layer=-1):
    nx, ny, _ = color_id.shape
    s, off = 1.0 - gap, gap / 2.0
    layer = color_id[:, :, z_layer]
    for x in range(nx):
        for y in range(ny):
            cid = int(layer[x, y])
            if cid < 0:
                continue
            ax.add_patch(Rectangle(
                (x + off, y + off), s, s,
                facecolor=PALETTE[cid % len(PALETTE)],
                edgecolor="k",
                linewidth=0.6,
            ))
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    apply_font(ax, title_pad=TITLE_PAD_2D)


def plot_cells(color_id, gap=0.18):
    fig3d = plt.figure()
    ax3d = fig3d.add_subplot(111, projection="3d")
    plot_3d(ax3d, color_id, gap=gap)
    fig3d.tight_layout()
    fig3d.savefig("task_split.png", dpi=200, bbox_inches="tight", pad_inches=0.4)

    fig2d = plt.figure()
    ax2d = fig2d.add_subplot(111)
    plot_top(ax2d, color_id, gap=gap, z_layer=-1)
    ax2d.set_title("Podział zadań - widok z góry")
    apply_font(ax2d, title_pad=TITLE_PAD_2D)
    fig2d.tight_layout()
    fig2d.savefig("task_split_top.png", dpi=200, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    plot_cells(color_id)