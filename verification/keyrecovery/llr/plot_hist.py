import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


def plot_histogram(csv_path: Path, pdf_path: Path, bins: int) -> None:
    data = pd.read_csv(csv_path, header=None, names=["bits"])
    values = data["bits"].to_numpy()

    stats = {
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
    }
    total = len(values)
    min_pct = float(np.count_nonzero(values == stats["min"]) * 100.0 / total)
    max_pct = float(np.count_nonzero(values == stats["max"]) * 100.0 / total)

    weights = np.full_like(values, 100.0 / len(values))
    legend_entries = [
        (f"min   : {stats['min']:.2f} ({min_pct:.2f}%)", "#222222"),
        (f"max   : {stats['max']:.2f} ({max_pct:.2f}%)", "#222222"),
        (f"mean  : {stats['mean']:.2f}", "#d62728"),
        (f"median: {stats['median']:.2f}", "#ff7f0e"),
    ]

    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.hist(
        values,
        bins=bins,
        weights=weights,
        color="#1f77b4",
        edgecolor="black",
    )
    mean_line = ax.axvline(stats["mean"], color="#d62728", linestyle="-", linewidth=1.4, label="mean")
    median_line = ax.axvline(
        stats["median"], color="#ff7f0e", linestyle="--", linewidth=1.4, label="median"
    )
    ax.set_yscale("log")
    ax.set_xlabel("Bits of information")
    ax.set_ylabel("Percentage (%) (log scale)")
    stats_handles = [
        Line2D([], [], linestyle="", linewidth=0, marker=None, label=label)
        for label, _ in legend_entries
    ]
    legend = ax.legend(
        handles=stats_handles,
        labels=[label for label, _ in legend_entries],
        loc="upper right",
        frameon=True,
        fancybox=True,
        handlelength=0.0,
        handletextpad=0.2,
    )
    for text, (_, color) in zip(legend.get_texts(), legend_entries):
        text.set_ha("left")
        text.set_color(color)
    legend._legend_box.align = "left"  # type: ignore[attr-defined]
    ax.grid(True, which="both", axis="both", linestyle="--", linewidth=0.6, alpha=0.6)
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10, subs=[1.0, 2.0, 5.0], numticks=20))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:.0e}" if y > 0 else "0"))
    ax.yaxis.set_minor_locator(ticker.NullLocator())

    fig.tight_layout()
    fig.savefig(pdf_path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot histogram from LLR experiment CSV output")
    parser.add_argument("csv", type=Path, help="Input CSV file containing one sample per line")
    parser.add_argument("pdf", nargs="?", type=Path, help="Output PDF path (defaults to CSV name with .pdf)")
    parser.add_argument("--bins", type=int, default=40, help="Number of histogram bins")
    args = parser.parse_args()

    output = args.pdf if args.pdf else args.csv.with_suffix(".pdf")
    plot_histogram(args.csv, output, args.bins)


if __name__ == "__main__":
    main()
