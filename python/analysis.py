import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =========================
# Загрузка основной статистики
# =========================

df = pd.read_csv("summary.csv")

# =========================
# ГРАФИКИ МАТЕМАТИЧЕСКИХ ОЖИДАНИЙ
# =========================

# ---------- mean_inside ----------
plt.figure(figsize=(8, 5))

for T in sorted(df["T"].unique()):
    part = df[df["T"] == T]
    plt.plot(
        part["Vmax"],
        part["mean_inside"],
        marker='o',
        label=f"T = {T}"
    )

plt.xlabel("Vmax")
plt.ylabel("Среднее время взаимодействия")
plt.title("Зависимость математического ожидания взаимодействия от Vmax")
plt.grid(True)
plt.legend()

plt.savefig("mean_inside_vs_vmax2d.png")
plt.show()

# ---------- mean_outside ----------
plt.figure(figsize=(8, 5))

for T in sorted(df["T"].unique()):
    part = df[df["T"] == T]
    plt.plot(
        part["Vmax"],
        part["mean_outside"],
        marker='o',
        label=f"T = {T}"
    )

plt.xlabel("Vmax")
plt.ylabel("Среднее время расхождения")
plt.title("Зависимость математического ожидания расхождения от Vmax")
plt.grid(True)
plt.legend()

plt.savefig("mean_outside_vs_vmax2d.png")
plt.show()


# =========================
# ФУНКЦИИ РАСПРЕДЕЛЕНИЯ ВЕРОЯТНОСТЕЙ
# =========================



def plot_cdf(filename, title, outname):

    data = pd.read_csv(filename, header=None)[0]

    data = np.sort(data)

    y = np.arange(1, len(data) + 1) / len(data)

    plt.figure(figsize=(8, 5))

    plt.plot(data, y)

    plt.xlabel("τ")
    plt.ylabel("F(τ)")
    plt.title(title)

    plt.grid(True)

    plt.savefig(outname)

    plt.show()



plot_cdf(
    "inside_Vmax_5.0_T_5.0.csv",
    "Функция распределения взаимодействия",
    "cdf_inside55.png"
)

plot_cdf(
    "outside_Vmax_5.0_T_5.0.csv",
    "Функция распределения расхождения",
    "cdf_outside55.png"
)


df2d = pd.read_csv("summary.csv")
df3d = pd.read_csv("summary3D.csv")


g2_in = df2d.groupby("Vmax")["count_inside"].mean()
g3_in = df3d.groupby("Vmax")["count_inside"].mean()

g2_out = df2d.groupby("Vmax")["count_outside"].mean()
g3_out = df3d.groupby("Vmax")["count_outside"].mean()

# =========================
# 1. COUNT INSIDE
# =========================
plt.figure(figsize=(8,5))
plt.plot(g2_in.index, g2_in.values, marker='o', label="2D")
plt.plot(g3_in.index, g3_in.values, marker='o', label="3D")

plt.title("Количество взаимодействий (inside)")
plt.xlabel("Vmax")
plt.ylabel("Count inside")
plt.grid(True)
plt.legend()
plt.savefig("count_inside_2d_3d.png")
plt.show()

# =========================
# 2. COUNT OUTSIDE
# =========================
plt.figure(figsize=(8,5))
plt.plot(g2_out.index, g2_out.values, marker='o', label="2D")
plt.plot(g3_out.index, g3_out.values, marker='o', label="3D")

plt.title("Количество расхождений (outside)")
plt.xlabel("Vmax")
plt.ylabel("Count outside")
plt.grid(True)
plt.legend()
plt.savefig("count_outside_2d_3d.png")
plt.show()