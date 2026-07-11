#!/usr/bin/env python3
"""Reusable: render a static image from a self-contained plotly `write_html` file by
extracting the figure JSON from its `Plotly.newPlot(uuid, [data], {layout}, {config})`
call and handing it to kaleido. Headless, high-res, respects the script's camera.

Usage: plotly_html_to_static.py IN.html OUT.png [--scale S] [--drop-title]
"""
import argparse
import json
import sys

import plotly.io as pio


def _scan_balanced(s, start, open_ch, close_ch):
    """Index of the close bracket matching s[start]==open_ch, string/escape aware."""
    depth = 0
    in_str = esc = False
    for i in range(start, len(s)):
        c = s[i]
        if in_str:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == '"':
                in_str = False
        elif c == '"':
            in_str = True
        elif c == open_ch:
            depth += 1
        elif c == close_ch:
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("unbalanced brackets")


def extract_figure(html_path):
    text = open(html_path, encoding="utf-8").read()
    idx = text.rfind("Plotly.newPlot(")            # the real (data-bearing) call is last
    if idx < 0:
        raise ValueError("no Plotly.newPlot call found")
    sub = text[idx:]
    lb = sub.index("[")                            # start of data array (after the uuid)
    de = _scan_balanced(sub, lb, "[", "]")
    data = json.loads(sub[lb:de + 1])
    ob = sub.index("{", de + 1)                    # start of layout object
    le = _scan_balanced(sub, ob, "{", "}")
    layout = json.loads(sub[ob:le + 1])
    return pio.from_json(json.dumps({"data": data, "layout": layout}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("html")
    ap.add_argument("out")
    ap.add_argument("--scale", type=float, default=2.0)
    ap.add_argument("--width", type=int, default=None)
    ap.add_argument("--height", type=int, default=None)
    ap.add_argument("--drop-title", action="store_true")
    ap.add_argument("--drop-text-traces", action="store_true",
                    help="remove Scatter3d traces in text-only mode (declutter labels)")
    ap.add_argument("--zoom", type=float, default=None,
                    help="scale the scene camera eye vector (<1 zooms in)")
    ap.add_argument("--colorbar-title", default=None,
                    help="override any trace colorbar title (e.g. method-neutral label)")
    a = ap.parse_args()

    fig = extract_figure(a.html)
    if a.drop_title:
        fig.update_layout(title=None, margin=dict(l=0, r=0, t=0, b=0))
    if a.drop_text_traces:
        fig.data = tuple(t for t in fig.data
                         if not (getattr(t, "mode", None) == "text"))
    if a.zoom is not None:
        cam = (fig.layout.scene.camera or {})
        eye = cam.eye if cam and cam.eye else None
        if eye and eye.x is not None:
            fig.layout.scene.camera.eye = dict(
                x=eye.x * a.zoom, y=eye.y * a.zoom, z=eye.z * a.zoom)
    if a.colorbar_title is not None:
        for t in fig.data:
            mk = getattr(t, "marker", None)
            if mk is not None and getattr(mk, "colorbar", None) is not None \
                    and mk.colorbar.title is not None \
                    and mk.colorbar.title.text is not None:
                mk.colorbar.title.text = a.colorbar_title
    kw = dict(scale=a.scale)
    if a.width:
        kw["width"] = a.width
    if a.height:
        kw["height"] = a.height
    fig.write_image(a.out, **kw)
    sz = __import__("os").path.getsize(a.out) / 1024
    print(f"wrote {a.out}  ({sz:.0f} KB)  traces={len(fig.data)}", file=sys.stderr)


if __name__ == "__main__":
    main()
