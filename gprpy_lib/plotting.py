import matplotlib.pyplot as plt
import numpy as np

def plot_profile(data, profilePos, twtt, velocity=None, maxTopo=None, minTopo=None,
                 color="gray", contrast=1.0, yrng=None, xrng=None, asp=None):
    """
    Plots the GPR profile data.

    INPUT:
    data         GPR data matrix
    profilePos   Array of profile positions
    twtt         Array of two-way travel times
    velocity     (Optional) Radar wave velocity in m/ns
    maxTopo      (Optional) Maximum topography value
    minTopo      (Optional) Minimum topography value
    color        Colormap for the plot (default: "gray")
    contrast     Contrast enhancement factor (default: 1.0)
    yrng         Y-axis range (default: None)
    xrng         X-axis range (default: None)
    asp          Aspect ratio (default: None)
    """
    dx = profilePos[3] - profilePos[2]
    dt = twtt[3] - twtt[2]
    stdcont = np.nanmax(np.abs(data)[:])

    fig, ax = plt.subplots()

    if velocity is None:
        im = ax.imshow(data, cmap=color,
                       extent=[min(profilePos) - dx / 2.0, max(profilePos) + dx / 2.0,
                               max(twtt) + dt / 2.0, min(twtt) - dt / 2.0],
                       aspect="auto", vmin=-stdcont / contrast, vmax=stdcont / contrast)
        ax.set_ylabel("time [ns]")
        ax.invert_yaxis()
        if yrng is not None:
            yrng = [np.max(yrng), np.min(yrng)]
        else:
            yrng = [np.max(twtt), np.min(twtt)]
    elif maxTopo is None:
        depth = twtt * velocity / 2.0
        dy = dt * velocity
        im = ax.imshow(data, cmap=color,
                       extent=[min(profilePos) - dx / 2.0, max(profilePos) + dx / 2.0,
                               max(depth) + dy / 2.0, min(depth) - dy / 2.0],
                       aspect="auto", vmin=-stdcont / contrast, vmax=stdcont / contrast)
        ax.set_ylabel("depth [m]")
        ax.invert_yaxis()
        if yrng is not None:
            yrng = [np.max(yrng), np.min(yrng)]
        else:
            yrng = [np.max(depth), np.min(depth)]
    else:
        depth = twtt * velocity / 2.0
        dy = dt * velocity
        im = ax.imshow(data, cmap=color,
                       extent=[min(profilePos) - dx / 2.0, max(profilePos) + dx / 2.0,
                               minTopo - max(depth) - dy / 2.0,
                               maxTopo - min(depth) + dy / 2.0],
                       aspect="auto", vmin=-stdcont / contrast, vmax=stdcont / contrast)
        ax.set_ylabel("elevation [m]")
        if yrng is None:
            yrng = [minTopo - np.max(depth), maxTopo - np.min(depth)]

    if xrng is None:
        xrng = [min(profilePos), max(profilePos)]

    if yrng is not None:
        ax.set_ylim(yrng)

    if xrng is not None:
        ax.set_xlim(xrng)

    if asp is not None:
        ax.set_aspect(asp)

    ax.set_xlabel("profile position [m]")
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    fig.colorbar(im, ax=ax, label='Amplitude')

    plt.show()
