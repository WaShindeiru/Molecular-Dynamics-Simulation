import re
from datetime import datetime
import plotly.graph_objects as go

LOG_PATH = "/mnt/7E442D59442D1585/md/optimization_attempt/threads_investigation/v3/perf.log"
KEEP_ITERS = {964}

if __name__ == "__main__":
    events = []

    with open(LOG_PATH) as f:
        current_iter = None
        for line in f:
            m = re.match(r'\[PERF (\S+)\] (.*)', line.strip())
            if not m:
                continue

            t = datetime.strptime(m.group(1), "%H:%M:%S.%f")
            label = m.group(2)

            iter_m = re.search(r'iter (\d+)', label)
            if iter_m:
                current_iter = int(iter_m.group(1))

            if current_iter in KEEP_ITERS:
                events.append((t, label))

    if not events:
        print("No matching events found.")
        exit(1)

    t0 = events[0][0]
    positions = [(t - t0).total_seconds() * 1000 for t, _ in events]
    labels = [label for _, label in events]

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=positions,
        y=[0] * len(positions),
        mode='markers',
        marker=dict(symbol='line-ns', size=16, line=dict(width=1)),
        hovertext=labels,
        hoverinfo='text+x',
    ))

    annotations = [
        dict(
            x=x, y=0,
            text=label,
            textangle=-90,
            showarrow=False,
            xanchor='left',
            yanchor='bottom',
            yshift=12,
        )
        for x, label in zip(positions, labels)
    ]

    fig.update_layout(
        xaxis_title="time (ms)",
        yaxis=dict(visible=False, range=[-1, 5]),
        title=f"log events — iters {sorted(KEEP_ITERS)}",
        height=600,
        annotations=annotations,
    )

    fig.show()
