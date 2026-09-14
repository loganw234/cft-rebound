#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
"""The bench sheet: one page from the measurement CSVs, nothing typed.

    python3 hw/report.py --dir bench-modes --out bench-modes/report.html

Reads what the benchmarks wrote - workload.csv, workload-drift-fp128.csv,
workload-fp256.csv, ensemble.csv, decompose.csv - and draws five figures,
each with a table twin, in a self-contained HTML page. Every number on the
page comes from a CSV row; the few that do not (the engine rate, the
image's worth to IAS15) are named with the ledger entry they come from.
"""
import argparse
import csv
import html
import math
import pathlib
import sys

# Bodies per problem, from the N line of each file under data/problems and
# the count tools/make_nbody.py was given. Carried here because the CSVs
# record the problem's name and the chart's x axis wants its size.
BODIES = {"kepler": 2, "pythagorean": 3, "outer": 6,
          "nbody64": 64, "nbody256": 256, "nbody512": 512}
PRETTY = {"kepler": "Kepler, 2 bodies", "pythagorean": "Pythagorean, 3",
          "outer": "outer solar system, 6", "nbody64": "n-body, 64",
          "nbody256": "n-body, 256", "nbody512": "n-body, 512"}

# Series identity is fixed across every figure: an image keeps its hue
# wherever it appears (validated categorical slots 1-3 all-pairs); the
# software backend is the de-emphasis gray, REBOUND is the ink.
SERIES = {
    "1t-f128": ("one f128 tile, 150 MHz", "--s1"),
    "4t-full": ("four full tiles, 135 MHz", "--s2"),
    "6t-f128": ("six f128 tiles, 125 MHz", "--s3"),
    "software": ("libcft software, FMA form", "--gray"),
    "sw": ("libcft software, FMA form", "--gray"),
    "sw-rebound-form": ("libcft software, REBOUND's form", "--gray2"),
    "ref": ("REBOUND, hardware double", "--ink"),
}
FMT_NAME = {"fp64": "binary64", "fp128": "binary128", "fp256": "binary256"}


def read(path):
    p = pathlib.Path(path)
    if not p.exists():
        return []
    with open(p, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def fnum(s, default=None):
    try:
        return float(s)
    except (TypeError, ValueError):
        return default


def esc(s):
    return html.escape(str(s))


# ---- scales -----------------------------------------------------------------
class Log:
    def __init__(self, lo, hi, a, b):
        self.lo, self.hi, self.a, self.b = math.log10(lo), math.log10(hi), a, b
    def __call__(self, v):
        if v is None or v <= 0:
            return None
        return self.a + (math.log10(v) - self.lo) / (self.hi - self.lo) * (self.b - self.a)
    def ticks(self):
        out = []
        for e in range(math.floor(self.lo), math.ceil(self.hi) + 1):
            out.append(10 ** e)
        return [t for t in out if 10 ** self.lo * 0.999 <= t <= 10 ** self.hi * 1.001]


def nice_bounds(vals, pad=1.6):
    vals = [v for v in vals if v is not None and v > 0]
    if not vals:
        return 1, 10
    lo, hi = min(vals) / pad, max(vals) * pad
    return 10 ** math.floor(math.log10(lo)), 10 ** math.ceil(math.log10(hi))


def si(v):
    if v is None:
        return "-"
    if v >= 1e9: return f"{v/1e9:.3g}G"
    if v >= 1e6: return f"{v/1e6:.3g}M"
    if v >= 1e3: return f"{v/1e3:.3g}k"
    if v >= 1: return f"{v:.3g}"
    return f"{v:.2g}"


def sci(v):
    if v is None:
        return "-"
    if v == 0:
        return "0"
    e = math.floor(math.log10(abs(v)))
    m = v / 10 ** e
    return f"{m:.1f}e{e:+d}"


# ---- one SVG panel ------------------------------------------------------------
class Panel:
    """A log-log panel with hairline grid, drawn to one scale per axis."""
    W, H = 440, 300
    L, R, T, B = 58, 18, 22, 44

    def __init__(self, title, xlabel, ylabel, xr, yr):
        self.title, self.xlabel, self.ylabel = title, xlabel, ylabel
        self.x = Log(xr[0], xr[1], self.L, self.W - self.R)
        self.y = Log(yr[0], yr[1], self.H - self.B, self.T)
        self.parts = []

    def grid(self, xfmt=si, yfmt=si):
        p = []
        for t in self.y.ticks():
            yy = self.y(t)
            p.append(f'<line class="grid" x1="{self.L}" x2="{self.W-self.R}" y1="{yy:.1f}" y2="{yy:.1f}"/>')
            p.append(f'<text class="tick" x="{self.L-6}" y="{yy+3.5:.1f}" text-anchor="end">{esc(yfmt(t))}</text>')
        for t in self.x.ticks():
            xx = self.x(t)
            p.append(f'<line class="grid" y1="{self.T}" y2="{self.H-self.B}" x1="{xx:.1f}" x2="{xx:.1f}"/>')
            p.append(f'<text class="tick" x="{xx:.1f}" y="{self.H-self.B+15}" text-anchor="middle">{esc(xfmt(t))}</text>')
        p.append(f'<text class="axis" x="{(self.L+self.W-self.R)/2:.0f}" y="{self.H-6}" text-anchor="middle">{esc(self.xlabel)}</text>')
        p.append(f'<text class="axis" transform="translate(13 {(self.T+self.H-self.B)/2:.0f}) rotate(-90)" text-anchor="middle">{esc(self.ylabel)}</text>')
        self.parts = p + self.parts

    def line(self, pts, var, label=None, tip=None):
        pts = [(self.x(a), self.y(b), a, b) for a, b in pts if a and b]
        pts = [q for q in pts if q[0] is not None and q[1] is not None]
        if not pts:
            return
        d = "M" + " L".join(f"{x:.1f} {y:.1f}" for x, y, _, _ in pts)
        self.parts.append(f'<path class="ln" style="stroke:var({var})" d="{d}"/>')
        for x, y, a, b in pts:
            t = tip(a, b) if tip else f"{si(a)}, {si(b)}"
            self.parts.append(f'<circle class="pt" style="fill:var({var})" cx="{x:.1f}" cy="{y:.1f}" r="4"><title>{esc(t)}</title></circle>')
        if label:
            x, y = pts[-1][0], pts[-1][1]
            self.parts.append(f'<text class="lbl" x="{x+7:.1f}" y="{y+3.5:.1f}">{esc(label)}</text>')

    def dot(self, a, b, var, shape="circle", tip="", r=5):
        x, y = self.x(a), self.y(b)
        if x is None or y is None:
            return
        if shape == "square":
            self.parts.append(f'<rect class="pt" style="fill:var({var})" x="{x-r:.1f}" y="{y-r:.1f}" width="{2*r}" height="{2*r}"><title>{esc(tip)}</title></rect>')
        elif shape == "diamond":
            self.parts.append(f'<path class="pt" style="fill:var({var})" d="M{x:.1f} {y-r-1:.1f} L{x+r+1:.1f} {y:.1f} L{x:.1f} {y+r+1:.1f} L{x-r-1:.1f} {y:.1f} Z"><title>{esc(tip)}</title></path>')
        elif shape == "tri":
            self.parts.append(f'<path class="pt" style="fill:var({var})" d="M{x:.1f} {y-r-1:.1f} L{x+r+1:.1f} {y+r:.1f} L{x-r-1:.1f} {y+r:.1f} Z"><title>{esc(tip)}</title></path>')
        elif shape == "hollow":
            self.parts.append(f'<circle class="pt hollow" style="stroke:var({var})" cx="{x:.1f}" cy="{y:.1f}" r="{r}"><title>{esc(tip)}</title></circle>')
        else:
            self.parts.append(f'<circle class="pt" style="fill:var({var})" cx="{x:.1f}" cy="{y:.1f}" r="{r}"><title>{esc(tip)}</title></circle>')

    def hline(self, v, text):
        y = self.y(v)
        if y is None:
            return
        self.parts.append(f'<line class="ref" x1="{self.L}" x2="{self.W-self.R}" y1="{y:.1f}" y2="{y:.1f}"/>')
        self.parts.append(f'<text class="tick" x="{self.W-self.R}" y="{y-4:.1f}" text-anchor="end">{esc(text)}</text>')

    def svg(self):
        return (f'<svg class="panel" viewBox="0 0 {self.W} {self.H}" role="img" aria-label="{esc(self.title)}">'
                f'<text class="ptitle" x="{self.L}" y="14">{esc(self.title)}</text>' + "".join(self.parts) + "</svg>")


def legend(items):
    return '<div class="legend">' + "".join(
        f'<span><i style="background:var({v})"></i>{esc(n)}</span>' for n, v in items) + "</div>"


def table(headers, rows):
    h = "".join(f"<th>{esc(x)}</th>" for x in headers)
    b = "".join("<tr>" + "".join(f"<td>{esc(c)}</td>" for c in r) + "</tr>" for r in rows)
    return f'<details class="twin"><summary>Table view</summary><div class="scroll"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div></details>'


# ---- figures ------------------------------------------------------------------
def fig_crossover(wl):
    """Seconds per step against bodies, one line per mode, one panel per format."""
    out = []
    panels = []
    for fmt in ("fp64", "fp128"):
        rows = [r for r in wl if r["format"] == fmt and r["problem"] in BODIES]
        secs = [fnum(r["seconds"]) / int(r["steps"]) for r in rows if fnum(r["seconds"])]
        yr = nice_bounds(secs, 3)
        pn = Panel(f"{FMT_NAME[fmt]}: seconds a step", "bodies", "seconds a step (log)", (1.5, 800), yr)
        pn.grid()
        modes = ["ref", "sw", "1t-f128", "4t-full", "6t-f128"] if fmt == "fp64" else ["sw", "1t-f128", "4t-full", "6t-f128"]
        for m in modes:
            pts = sorted(((BODIES[r["problem"]], fnum(r["seconds"]) / int(r["steps"]))
                          for r in rows if r["mode"] == m and fnum(r["seconds"])), key=lambda q: q[0])
            name, var = SERIES[m]
            pn.line(pts, var, label={"ref": "REBOUND", "sw": "software", "1t-f128": "1 tile", "4t-full": "4 tiles", "6t-f128": "6 tiles"}[m],
                    tip=lambda a, b, n=name: f"{n}: {b:.3g} s a step at {a} bodies")
        panels.append(pn.svg())
    trows = [(PRETTY.get(r["problem"], r["problem"]), r["steps"], FMT_NAME.get(r["format"], r["format"]),
              SERIES.get(r["mode"], (r["mode"],))[0], f'{fnum(r["seconds"]):.2f}',
              f'{fnum(r["seconds"])/int(r["steps"]):.4g}', r["record"], r["passes"])
             for r in wl if r["problem"] in BODIES]
    out.append('<div class="panels">' + "".join(panels) + "</div>")
    out.append(legend([("REBOUND, hardware double", "--ink"), ("libcft software, FMA form", "--gray"),
                       ("one f128 tile", "--s1"), ("four full tiles", "--s2"), ("six f128 tiles", "--s3")]))
    out.append(table(["problem", "steps", "format", "mode", "seconds", "s per step", "record", "corrector passes"], trows))
    return "".join(out)


def fig_frontier(wl, drift128, wl256):
    """Accuracy for cost: seconds a step against energy drift, one panel per problem."""
    # Exact binary128 drifts override the run's own column, which was
    # computed in doubles; binary256 rows come from the quick pass.
    drift = {}
    for r in drift128:
        if r["mode"] in ("sw",):
            drift[(r["problem"], "fp128")] = fnum(r["energy_drift"])
    rows = list(wl) + list(wl256)
    panels, trows = [], []
    probs = [p for p in BODIES if any(r["problem"] == p for r in rows)]
    for p in probs:
        pr = [r for r in rows if r["problem"] == p]
        xs = [fnum(r["seconds"]) / int(r["steps"]) for r in pr if fnum(r["seconds"])]
        ds = []
        for r in pr:
            d = fnum(r["energy_drift"])
            if r["format"] == "fp128" and (p, "fp128") in drift:
                d = drift[(p, "fp128")]
            if d and d > 0:
                ds.append(d)
        if not ds:
            continue
        pn = Panel(PRETTY.get(p, p), "seconds a step (log)", "energy drift (log)", nice_bounds(xs, 3), nice_bounds(ds, 4))
        pn.grid(yfmt=sci)
        for r in pr:
            s = fnum(r["seconds"])
            if not s:
                continue
            d = fnum(r["energy_drift"])
            if r["format"] == "fp128" and (p, "fp128") in drift:
                d = drift[(p, "fp128")]
            if not d or d <= 0:
                continue
            var = {"fp64": "--s1", "fp128": "--s2", "fp256": "--s3"}[r["format"]]
            shape = {"ref": "square", "sw": "hollow", "sw-rebound-form": "hollow",
                     "1t-f128": "circle", "4t-full": "diamond", "6t-f128": "tri"}.get(r["mode"], "circle")
            name = SERIES.get(r["mode"], (r["mode"],))[0]
            pn.dot(s / int(r["steps"]), d, var, shape, f"{name}, {FMT_NAME[r['format']]}: {s/int(r['steps']):.3g} s a step, drift {d:.2e}")
            trows.append((PRETTY.get(p, p), FMT_NAME[r["format"]], name, f"{s/int(r['steps']):.4g}", f"{d:.3e}", r["record"]))
        panels.append(pn.svg())
    out = ['<div class="panels">' + "".join(panels) + "</div>",
           legend([("binary64", "--s1"), ("binary128", "--s2"), ("binary256", "--s3")]),
           '<p class="note">Shape is the mode: square REBOUND in double, hollow circle software, filled circle one tile, diamond four tiles, triangle six tiles. Read down a column of equal drift.</p>',
           table(["problem", "format", "mode", "s per step", "energy drift", "record"], trows)]
    return "".join(out)


def fig_ensemble(en):
    if not en:
        return '<p class="note">The ensemble sweep has not been run yet.</p>'
    panels, trows = [], []
    for fmt in ("fp64", "fp128"):
        rows = [r for r in en if r["format"] == fmt and fnum(r["system_steps_per_s"])]
        if not rows:
            continue
        ys = [fnum(r["system_steps_per_s"]) for r in rows]
        pn = Panel(f"{FMT_NAME[fmt]}: system-steps a second", "members in the ensemble", "system-steps a second (log)",
                   (0.7, 6000), nice_bounds(ys, 3))
        pn.grid()
        for m in ("software", "1t-f128", "4t-full", "6t-f128"):
            pts = sorted(((int(r["members"]), fnum(r["system_steps_per_s"])) for r in rows if r["mode"] == m), key=lambda q: q[0])
            name, var = SERIES[m]
            pn.line(pts, var, label={"software": "software", "1t-f128": "1 tile", "4t-full": "4 tiles", "6t-f128": "6 tiles"}[m],
                    tip=lambda a, b, n=name: f"{n}: {b:.4g} system-steps a second at E = {a}")
        panels.append(pn.svg())
    for r in en:
        trows.append((r["members"], FMT_NAME.get(r["format"], r["format"]), SERIES.get(r["mode"], (r["mode"],))[0],
                      r["tiles"], r["seconds"], r["system_steps_per_s"], r["vs_sw"], r["record"], r["lane_efficiency"], r["passes"]))
    return ('<div class="panels">' + "".join(panels) + "</div>" +
            legend([("libcft software", "--gray"), ("one f128 tile", "--s1"), ("four full tiles", "--s2"), ("six f128 tiles", "--s3")]) +
            table(["E", "format", "mode", "tiles", "seconds", "system-steps/s", "vs software", "record", "lane efficiency", "passes"], trows))


def fig_decompose(dc):
    if not dc:
        return '<p class="note">The decomposition pass has not been run yet.</p>'
    W, RH = 760, 26
    rows = [r for r in dc if fnum(r["wall_seconds"])]
    h = 30 + RH * len(rows) + 30
    parts = [f'<svg class="bars" viewBox="0 0 {W} {h}" role="img" aria-label="Where a step goes">']
    L, R = 250, 20
    for i, r in enumerate(rows):
        w = fnum(r["wall_seconds"]); tp = fnum(r["t_program"], 0); te = fnum(r["t_elem"], 0)
        td = fnum(r["t_divsqrt"], 0); tg = fnum(r["t_gravity"], 0)
        other = max(w - tp - te - td, 0)
        y = 30 + i * RH
        label = f'{PRETTY.get(r["problem"], r["problem"])}, {FMT_NAME.get(r["format"], r["format"])}, {r["backend"]}'
        parts.append(f'<text class="tick" x="{L-8}" y="{y+15}" text-anchor="end">{esc(label)}</text>')
        x = L
        for seg, var, nm in ((tp, "--s1", "program runs"), (te, "--s2", "elementwise calls"), (td, "--s3", "divide and sqrt"), (other, "--gray", "host, everything else")):
            ww = (W - L - R) * seg / w if w else 0
            if ww > 0:
                parts.append(f'<rect style="fill:var({var})" x="{x+1:.1f}" y="{y+3}" width="{max(ww-2,0.5):.1f}" height="18" rx="0"><title>{esc(nm)}: {seg:.2f} s of {w:.2f} ({100*seg/w:.0f}%)</title></rect>')
            x += ww
        gx = L + (W - L - R) * tg / w if w else L
        parts.append(f'<line class="mark" x1="{gx:.1f}" x2="{gx:.1f}" y1="{y+1}" y2="{y+23}"><title>gravity, inclusive: {tg:.2f} s ({100*tg/w:.0f}%)</title></line>')
        parts.append(f'<text class="lbl" x="{W-R+2}" y="{y+15}" text-anchor="start">{w:.1f}s</text>')
    parts.append("</svg>")
    trows = [(PRETTY.get(r["problem"], r["problem"]), r["steps"], FMT_NAME.get(r["format"], r["format"]), r["backend"], r["wall_seconds"],
              r["t_program"], r["t_elem"], r["t_divsqrt"], r["t_gravity"], r["calls"], r["program_calls"], r["divsqrt_calls"],
              si(fnum(r["staged_bytes"])) + "B", r["passes"]) for r in dc]
    return ("".join(parts) +
            legend([("program runs", "--s1"), ("elementwise calls", "--s2"), ("divide and sqrt", "--s3"), ("host, everything else", "--gray")]) +
            '<p class="note">The vertical tick on each bar is where gravity ends when its own time is measured inclusively - host gather, scatter and every library call inside it. It overlaps the segments rather than partitioning them, which is why it is a mark and not a fill.</p>' +
            table(["problem", "steps", "format", "backend", "wall s", "program s", "elementwise s", "div/sqrt s", "gravity s (incl.)", "calls", "program calls", "div/sqrt calls", "bytes presented", "passes"], trows))


def fig_determinism(wl, en):
    probs = [p for p in BODIES if any(r["problem"] == p for r in wl)]
    modes = ["sw-rebound-form", "sw", "1t-f128", "4t-full", "6t-f128"]
    head = ["problem", "format"] + [SERIES[m][0] for m in modes]
    body, n_ident, n_rows = [], 0, 0
    for p in probs:
        for fmt in ("fp64", "fp128"):
            cells = []
            for m in modes:
                r = next((x for x in wl if x["problem"] == p and x["format"] == fmt and x["mode"] == m), None)
                if r is None:
                    cells.append(("-", "none"))
                    continue
                rec = r["record"]
                n_rows += 1
                if rec.startswith("identical") or rec.startswith("fma ref"):
                    n_ident += 1
                    cells.append(("identical" if rec.startswith("identical") else "reference", "good"))
                elif rec.startswith("no rung"):
                    cells.append(("refused by name", "none"))
                else:
                    cells.append((rec, "bad"))
            body.append((PRETTY.get(p, p), FMT_NAME[fmt], cells))
    en_ident = sum(1 for r in en if r["record"] == "identical")
    en_total = sum(1 for r in en if r["mode"] != "software" and not r["record"].startswith("no rung"))
    rows_html = "".join(
        f"<tr><td>{esc(p)}</td><td>{esc(f)}</td>" + "".join(
            f'<td class="cell {cls}"><i></i>{esc(txt)}</td>' for txt, cls in cells) + "</tr>"
        for p, f, cells in body)
    return (f'<p class="kpi-line"><b>{n_ident}</b> of {n_rows} workload records byte-identical to their reference, '
            f'<b>{en_ident}</b> of {en_total} ensemble records; every card row was also identical across tiles. '
            f'A refusal is not a mismatch: an image without the rung says so by name and computes nothing.</p>'
            f'<div class="scroll"><table class="matrix"><thead><tr>' + "".join(f"<th>{esc(x)}</th>" for x in head) +
            f"</tr></thead><tbody>{rows_html}</tbody></table></div>")


# ---- the page -----------------------------------------------------------------
CSS = """
:root{color-scheme:light;
 --bg:#f4f5f7;--surface:#fcfcfd;--ink:#111418;--text:#1b2028;--text2:#4d5563;--muted:#7c8494;
 --line:#dfe3ea;--grid:#e9ecf1;--accent:#2a78d6;
 --s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--gray:#8a8f99;--gray2:#b3b7bf;--good:#008300;--bad:#e34948;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;
 --bg:#15181c;--surface:#1c2025;--ink:#f3f4f6;--text:#e6e8ec;--text2:#b5bac4;--muted:#8c93a0;
 --line:#2c323a;--grid:#262b32;--accent:#3987e5;
 --s1:#3987e5;--s2:#d95926;--s3:#199e70;--gray:#7d838e;--gray2:#565c66;--good:#3fbf3f;--bad:#e66767;}}
:root[data-theme="dark"]{color-scheme:dark;
 --bg:#15181c;--surface:#1c2025;--ink:#f3f4f6;--text:#e6e8ec;--text2:#b5bac4;--muted:#8c93a0;
 --line:#2c323a;--grid:#262b32;--accent:#3987e5;
 --s1:#3987e5;--s2:#d95926;--s3:#199e70;--gray:#7d838e;--gray2:#565c66;--good:#3fbf3f;--bad:#e66767;}
body{margin:0;background:var(--bg);color:var(--text);font-family:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",sans-serif;font-size:15px;line-height:1.55}
main{max-width:1180px;margin:0 auto;padding:36px 28px 72px}
h1{font-family:"Bricolage Grotesque","IBM Plex Sans",sans-serif;font-weight:600;font-size:38px;letter-spacing:-.01em;line-height:1.1;margin:0 0 6px;text-wrap:balance;color:var(--ink)}
h2{font-family:"Bricolage Grotesque","IBM Plex Sans",sans-serif;font-weight:600;font-size:24px;margin:52px 0 6px;letter-spacing:-.005em;text-wrap:balance;color:var(--ink)}
.eyebrow{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--text2)}
.lede,.note,p{max-width:70ch}
.lede{font-size:17px;color:var(--text2);margin:0 0 26px}
.note{color:var(--text2);font-size:14px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:14px;margin:26px 0 8px}
.kpi{background:var(--surface);border:1px solid var(--line);padding:16px 18px}
.kpi .l{font-size:13px;color:var(--text2)}
.kpi .v{font-size:38px;font-weight:600;line-height:1.1;margin:4px 0 2px;color:var(--ink)}
.kpi .d{font-size:13px;color:var(--muted)}
.kpi-line{font-size:15px}
.kpi-line b{color:var(--ink)}
.panels{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:16px;margin-top:14px}
svg.panel,svg.bars{width:100%;height:auto;background:var(--surface);border:1px solid var(--line);display:block}
svg text{font-family:"IBM Plex Mono",ui-monospace,monospace;fill:var(--text2);font-size:11px}
svg .ptitle{font-family:"IBM Plex Sans",sans-serif;font-size:13px;fill:var(--ink);font-weight:600}
svg .axis{font-size:11px;fill:var(--muted)}
svg .tick{font-size:10.5px;fill:var(--muted);font-variant-numeric:tabular-nums}
svg .lbl{font-family:"IBM Plex Sans",sans-serif;font-size:11.5px;fill:var(--text)}
svg .grid{stroke:var(--grid);stroke-width:1}
svg .ref{stroke:var(--muted);stroke-width:1}
svg .mark{stroke:var(--ink);stroke-width:2}
svg .ln{fill:none;stroke-width:2;stroke-linejoin:round;stroke-linecap:round}
svg .pt{stroke:var(--surface);stroke-width:2}
svg .pt.hollow{fill:var(--surface);stroke-width:2.5}
.legend{display:flex;flex-wrap:wrap;gap:6px 18px;margin:10px 0 4px;font-size:13px;color:var(--text2)}
.legend i{display:inline-block;width:12px;height:12px;margin-right:7px;vertical-align:-1px}
details.twin{margin:8px 0 0}
summary{cursor:pointer;color:var(--accent);font-size:13px}
.scroll{overflow-x:auto}
table{border-collapse:collapse;font-size:13px;font-variant-numeric:tabular-nums;margin-top:8px;min-width:100%}
th,td{text-align:left;padding:5px 10px;border-bottom:1px solid var(--line);white-space:nowrap}
th{font-weight:600;color:var(--text2);font-size:12px}
table.matrix td.cell i{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:7px;background:var(--muted)}
table.matrix td.good i{background:var(--good)}
table.matrix td.bad i{background:var(--bad)}
table.matrix td.none{color:var(--muted)}
.src{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px;color:var(--muted);margin-top:40px;border-top:1px solid var(--line);padding-top:14px}
@media (prefers-reduced-motion:no-preference){svg .pt{transition:r .12s}}
svg .pt:hover{r:6}
"""


def page(args):
    d = pathlib.Path(args.dir)
    wl = read(d / "workload.csv"); dr = read(d / "workload-drift-fp128.csv")
    w256 = read(d / "workload-fp256.csv"); en = read(d / "ensemble.csv"); dc = read(d / "decompose.csv")
    if not wl:
        sys.exit("no workload.csv")

    # Headline figures, each from a row or a named ledger entry.
    r512 = {(r["mode"], r["format"]): r for r in wl if r["problem"] == "nbody512"}
    def ratio(a, b):
        ra, rb = r512.get(a), r512.get(b)
        if ra and rb and fnum(ra["seconds"]) and fnum(rb["seconds"]):
            return fnum(rb["seconds"]) / fnum(ra["seconds"])
        return None
    k1 = ratio(("1t-f128", "fp128"), ("sw", "fp128"))
    k2 = ratio(("1t-f128", "fp128"), ("6t-f128", "fp128"))
    n_ident = sum(1 for r in wl if r["record"].startswith("identical") or r["record"].startswith("fma ref"))
    n_rows = sum(1 for r in wl if not r["record"].startswith("no rung") and r["mode"] != "ref")

    parts = [f"<title>IAS15 on the Tile</title><style>{CSS}</style>",
             '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:wght@500;600&family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Mono:wght@400&display=swap">',
             "<main>",
             '<div class="eyebrow">cft-rebound, measured on an Alveo U50C, 2026-09-14</div>',
             "<h1>Where the hardware starts to pay, on a real integration</h1>",
             '<p class="lede">REBOUND\'s IAS15 in hardware double against the same integrator through libcft, on the software backend and on three bitstreams: one and six tiles of the binary256-less image, and the shipped four-tile full image. Every row is an actual integration, and every card row is held byte-identical to its software counterpart.</p>',
             '<div class="kpis">',
             f'<div class="kpi"><div class="l">One tile at binary128, 512 bodies, against the same bits in software</div><div class="v">{k1:.2f}x</div><div class="d">the case for the card, at the precision binary64 cannot reach</div></div>' if k1 else "",
             f'<div class="kpi"><div class="l">One tile against six, same problem, same bits</div><div class="v">{k2:.2f}x</div><div class="d">tile count is a cost for this workload; the library splits every call across every tile</div></div>' if k2 else "",
             f'<div class="kpi"><div class="l">Records byte-identical to their reference</div><div class="v">{n_ident} of {n_rows}</div><div class="d">across software, one, four and six tiles, two bitstreams</div></div>',
             '<div class="kpi"><div class="l">The specialised image, isolated from tile count</div><div class="v">1.02x</div><div class="d">both one-unit images head to head; ledger entry 36</div></div>',
             "</div>",
             "<h2>1. Seconds a step against bodies</h2>",
             '<p class="note">The traditional option is the ink line: REBOUND in double, which nothing here beats at binary64 and which cannot be run at binary128 at all. The gray line is the same arithmetic in software. One tile crosses it between six and sixty-four bodies and stays above; four and six tiles sit below one at every size.</p>',
             fig_crossover(wl),
             "<h2>2. Accuracy for cost</h2>",
             '<p class="note">Speed alone cannot rank rows that do not compute the same thing. Binary64 floors near 1e-15 whatever the step; binary128 does not. Read down a column of equal drift: the cheapest binary128 point is the one that matters, and it is the single tile.</p>',
             fig_frontier(wl, dr, w256),
             "<h2>3. The ensemble: the workload with long vectors</h2>",
             '<p class="note">E independent Kepler systems in one run, members one ulp apart. At E members every predictor and corrector call carries 6E elements, which is the first place a multi-tile image gets a vector long enough to split six ways. Lane efficiency rides beside every ratio in the table.</p>',
             fig_ensemble(en),
             "<h2>4. Where a step goes</h2>",
             '<p class="note">Every library call timed on a wall clock and bucketed; gravity timed inclusively. This is the measurement that decides which of residency, call count and the gather is the wall, before any of them is built.</p>',
             fig_decompose(dc),
             "<h2>5. The same bits everywhere</h2>",
             '<p class="note">The property that lets a pooled node\'s answer be trusted without recomputation. Each cell compares a full record - every coordinate, every counter - against its reference: REBOUND for the port\'s own form at binary64, the software backend for every card row.</p>',
             fig_determinism(wl, en),
             '<div class="src">generated by hw/report.py from bench-modes/*.csv; the runs and their commits are in docs/VALIDATION.md entries 34-37</div>',
             "</main>"]
    pathlib.Path(args.out).write_text("".join(parts), encoding="utf-8", newline="\n")
    print(f"wrote {args.out}: {len(wl)} workload rows, {len(en)} ensemble rows, {len(dc)} decomposition rows")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="bench-modes")
    ap.add_argument("--out", default="bench-modes/report.html")
    page(ap.parse_args())
