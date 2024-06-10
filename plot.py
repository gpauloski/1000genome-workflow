from __future__ import annotations

import argparse
import pathlib
import sys
from typing import Sequence

import plotly.express as px
import plotly.figure_factory as ff
import polars
from plotly.subplots import make_subplots

STAGES = {
    "processing_chrom_parts": "Stage 1",
    "individuals": "Stage 1",
    "processing_2": "Stage 2",
    "merge": "Stage 2",
    "sifting": "Stage 3",
    "mutation_overlap": "Stage 4",
    "frequency": "Stage 5",
}

RUNS = {
    "gc-trial": "Baseline",
    "proxy-bench-trial": "ProxyFutures",
    "process-pool": "Baseline (ProcessPool)",
    "process-pool-proxystore": "ProxyFutures (ProcessPool)",
    "globus-compute": "Baseline",
    "globus-compute-proxystore": "ProxyFutures",
}


def load_data(filepaths: list[str]) -> polars.DataFrame:
    dfs = [polars.read_csv(filepath) for filepath in filepaths]
    return polars.concat(dfs)


def merge_stages(data: polars.DataFrame) -> polars.DataFrame:
    out = []
    for run in RUNS:
        run_data = data.filter(polars.col("name") == run)
        if len(run_data) == 0:
            continue
        run_start = run_data["start"].min()
        for stage in STAGES:
            stage_data = run_data.filter(polars.col("task") == stage)
            if len(stage_data) == 0:
                continue
            for row in stage_data.iter_rows(named=True):
                start = row["start"] - run_start
                end = row["end"] - run_start
                out.append(
                    {
                        "run": RUNS[run],
                        "Start": start,
                        "Finish": end,
                        "Task": STAGES[stage],
                    },
                )

    return polars.DataFrame(out)


def make_gantt(data, layout=None):
    runs = set(data.select(polars.col("run")).to_series().to_list())
    fig = make_subplots(
        cols=1,
        rows=len(runs),
        vertical_spacing=0.05,
        shared_xaxes=True,
    )

    max_x_range = 0
    colors = [0, 4, 1, 2, 3]
    for row, run in enumerate(runs):
        run_data = data.filter(polars.col("run") == run)
        max_x_range = max(max_x_range, run_data["Finish"].max())
        timeline_fig = ff.create_gantt(
            run_data.to_dicts(),
            index_col="Task",
            show_colorbar=False,
            group_tasks=True,
            bar_width=0.475,
            colors=[px.colors.qualitative.Bold[i] for i in colors],
        )
        for stage in range(5):
            x = run_data.filter(polars.col("Task") == f"Stage {stage+1}")["Start"].min()
            fig.add_annotation(
                x=x + (20 if stage != 3 else 30),
                y=5 - (stage + 1),
                align="left",
                text=stage + 1,
                font=dict(size=10, color="white"),
                showarrow=False,
                row=row + 1,
                col=1,
            )
        # Crazy hacks to add traces in right order for legend
        traces = list(timeline_fig.data)
        for stage in [f"Stage {i}" for i in range(5)]:
            for trace in traces:
                if trace.name == stage:
                    fig.add_trace(trace, row=row + 1, col=1)
                    traces.remove(trace)
                    break
        for trace in traces:
            fig.add_trace(trace, row=row + 1, col=1)
        fig.update_yaxes(
            title_text=run,
            showticklabels=False,
            ticks="",
            row=row + 1,
            col=1,
            # autorange='reversed',
        )
        if row + 1 == len(runs):
            fig.update_xaxes(
                title_text="Runtime (s)",
                range=[0, max_x_range + 0.1],
                row=row + 1,
                col=1,
            )

    fig.update_traces(showlegend=True)
    names = set()
    fig.for_each_trace(
        lambda trace: trace.update(showlegend=False)
        if trace.name == "" or trace.name in names
        else names.add(trace.name),
    )

    layout = {} if layout is None else layout
    fig.update_layout(
        legend=dict(
            title=None,
            yanchor="top",
            y=0.85,
            xanchor="left",
            x=1.01,
            tracegroupgap=1,
        ),
        **layout,
    )

    return fig


def get_deltas(data: polars.DataFrame) -> polars.DataFrame:
    processed = data.group_by("run", "Task").agg(
        polars.col("Start").min(),
        polars.col("Finish").max(),
    )
    processed = processed.select(
        polars.col("run").alias("Run"),
        polars.col("Task").alias("Stage"),
        (polars.col("Finish") - polars.col("Start")).alias("Runtime"),
    )
    for run in set(data.select(polars.col("run")).to_series().to_list()):
        run_data = data.filter(polars.col("run") == run)
        runtime = run_data["Finish"].max() - run_data["Start"].min()
        print(f"{run} runtime: {runtime:.2f} s")
    processed = processed.pivot(
        index="Stage",
        columns="Run",
        values=None,
    ).sort(polars.col("Stage"))
    
    baseline_name, proxystore_name = sorted(set(data['run']))
    
    processed = processed.with_columns(
        (polars.col(baseline_name) - polars.col(proxystore_name)).alias(
            "Delta (s)",
        ),
        (
            100
            * (polars.col(baseline_name) - polars.col(proxystore_name))
            / polars.col(baseline_name)
        ).alias("Delta (%)"),
    )
    return processed


def main(argv: Sequence[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]

    parser = argparse.ArgumentParser(description="Create a GANTT chart from run traces")
    parser.add_argument(
        "--baseline", required=True, help="path to output CSV file from baseline run",
    )
    parser.add_argument(
        "--proxystore",
        required=True,
        help="path to output CSV file from ProxyStore run",
    )
    parser.add_argument("--output", required=True, help="path to save GANTT chart to")
    args = parser.parse_args()

    raw_data = load_data([args.baseline, args.proxystore])
    data = merge_stages(raw_data)
    print(get_deltas(data))

    fig = make_gantt(data)

    fig_path = pathlib.Path(args.output)
    fig_path.parent.mkdir(exist_ok=True)
    fig.write_image(fig_path)

    print(f"Saved GANTT chart to {fig_path}")


if __name__ == "__main__":
    raise SystemExit(main())
