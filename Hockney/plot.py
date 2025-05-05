import matplotlib.pyplot as plt
import numpy as np


def plot_alternating_charges_polygon(sides, radius=1.0, ax=None):
    angles = np.linspace(0, 2 * np.pi, sides, endpoint=False)
    x = radius * np.cos(angles)
    y = radius * np.sin(angles)
    charges = np.array(["+" if i % 2 == 0 else "-" for i in range(sides)])

    if ax is None:
        fig, ax = plt.subplots()

    ax.set_title(f"{sides}-gon")
    for i in range(sides):
        color = "red" if charges[i] == "+" else "blue"
        ax.plot(x[i], y[i], "o", color=color)
        ax.text(
            x[i] * 1.1, y[i] * 1.1, charges[i], ha="center", va="center", fontsize=12
        )

    ax.plot(np.append(x, x[0]), np.append(y, y[0]), "k--", linewidth=0.5)  # outline
    ax.set_aspect("equal")
    ax.axis("off")


fig, axs = plt.subplots(1, 4, figsize=(16, 4))

sides_list = [4, 6, 8, 12]
for ax, sides in zip(axs, sides_list):
    plot_alternating_charges_polygon(sides, ax=ax)

plt.tight_layout()
plt.show()
