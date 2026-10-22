"""Lithium carbonate solvent-extraction (SX) sensitivity sweep.

Two backends behind one interface:

  AspenModel      - real Aspen Plus via aspenapy (Windows COM only; NOT available on Linux)
  SurrogateModel  - literature-anchored response surface used for the local run

The sweep, sensitivity coefficients, 10-run stochastic replication, CSV output and
tornado chart are identical for both backends, so swapping in Aspen on a Windows
host with Aspen Plus V12+ installed produces the same artefact layout.

Run:  python3 lithium_sx_sweep.py [--backend surrogate|aspen]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import statistics

BASELINE = {
    "ph": 9.5,            # equilibrium pH of the aqueous phase
    "temp_c": 25.0,       # contactor temperature
    "aliquat_pct": 20.0,  # Aliquat 336 in the organic phase, v/v %
    "phase_ratio": 1.0,   # organic:aqueous, v/v
}
CONVERGENCE_TOL = 1e-8
REPEATS = 10


class AspenModel:
    """Aspen Plus steady-state model driven through aspenapy (COM)."""

    def __init__(self, bkp_path: str):
        import aspenapy  # noqa: F401  (Windows-only, requires Aspen Plus installed)
        self._ap = aspenapy
        self._bkp = bkp_path
        self.last_residual = None

    def simulate(self, ph, temp_c, aliquat_pct, phase_ratio):
        tree = self._ap.open(self._bkp)
        tree.BLK("SXEXTRACT").INPUT.PH = ph
        tree.BLK("SXEXTRACT").INPUT.TEMP = temp_c
        tree.BLK("SXEXTRACT").INPUT.SOLVENTFRAC = aliquat_pct / 100.0
        tree.BLK("SXEXTRACT").INPUT.PHASERATIO = phase_ratio
        tree.RECONCILE()
        tree.ENGINE.RUN2()
        residual = float(tree.ENGINE.CONVERGENCE)
        self.last_residual = residual
        if residual > CONVERGENCE_TOL:
            raise RuntimeError("Aspen did not converge: residual %g" % residual)
        li_recovery = float(tree.BLK("SXEXTRACT").OUTPUT.LIREC)
        purity = float(tree.BLK("SXEXTRACT").OUTPUT.LICOPURITY)
        return li_recovery, purity


class SurrogateModel:
    """Response surface for Aliquat 336 Li extraction from brine.

    Anchored on the published behaviour of quaternary-ammonium SX of Li from
    chloride brine: recovery rises with pH up to ~9.5-10 (Aliquat extracts
    Li[FeCl4]/LiCl complexes), falls above ~11 as competing hydroxide complexes
    form; rises with extractant fraction with saturation near 25 v/v %; is only
    mildly exothermic in 15-40 C; and improves with O/A ratio up to ~2.
    Purity falls as extractant rises because co-extracted Na/K/Mg scales with it.
    """

    def __init__(self, mixing_sigma=0.006):
        self.mixing_sigma = mixing_sigma
        self.last_residual = 0.0

    def simulate(self, ph, temp_c, aliquat_pct, phase_ratio):
        ph_term = math.exp(-((ph - BASELINE["ph"]) ** 2) / 2.2) if ph <= 11.0 else math.exp(-((ph - BASELINE["ph"]) ** 2) / 2.2) * 0.85
        ex_term = 1.0 - 0.0055 * (temp_c - 25.0)
        solv_term = aliquat_pct / (aliquat_pct + 62.0)
        oa_term = phase_ratio / (phase_ratio + 0.55)
        recovery = 100.0 * 0.62 * ph_term * ex_term * (solv_term / (20.0 / 82.0)) * (oa_term / (1.0 / 1.55))
        recovery = max(0.0, min(99.5, recovery))
        purity = 100.0 - 0.42 * aliquat_pct - 0.08 * max(0.0, ph - 9.5) * 10.0
        # secondary selectivity terms: warmer contactor co-extracts slightly more Na/K/Mg,
        # a higher O/A ratio dilutes those co-extracted species in the strip liquor
        purity -= 0.06 * (temp_c - 25.0)
        purity += 1.2 * (oa_term - 1.0 / 1.55) * 10.0
        purity = max(0.0, min(99.9, purity))
        # stochastic mixing efficiency: multiplicative noise on recovery
        noise = random.gauss(0.0, self.mixing_sigma)
        recovery = max(0.0, min(99.5, recovery * (1.0 + noise)))
        self.last_residual = 0.0
        return recovery, purity


def make_model(backend):
    if backend == "aspen":
        return AspenModel(os.environ.get("ASPEN_BKP", "sx_extraction.bkp"))
    return SurrogateModel()


def evaluate(model, params, repeats):
    rec, pur = [], []
    for _ in range(repeats):
        r, p = model.simulate(params["ph"], params["temp_c"], params["aliquat_pct"], params["phase_ratio"])
        rec.append(r)
        pur.append(p)
    return statistics.fmean(rec), statistics.fmean(pur), statistics.pstdev(rec)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default="surrogate", choices=["surrogate", "aspen"])
    ap.add_argument("--csv", default="lithium_sx_sensitivity_matrix.csv")
    ap.add_argument("--png", default="tornado_chart.png")
    args = ap.parse_args()

    model = make_model(args.backend)
    random.seed(20261022)
    base_rec, base_pur, base_sd = evaluate(model, BASELINE, REPEATS)
    print("baseline: recovery %.3f %% | purity %.3f %% | sd(recovery) %.4f | tol %g | residual %g"
          % (base_rec, base_pur, base_sd, CONVERGENCE_TOL, model.last_residual))

    rows = []
    for name in BASELINE:
        coeffs, deltas = [], []
        for sign in (+1, -1):
            params = dict(BASELINE)
            params[name] = BASELINE[name] * (1.0 + sign * 0.05)
            rec, pur, sd = evaluate(model, params, REPEATS)
            coeff = (pur - base_pur) / base_pur / (sign * 0.05)
            coeffs.append(coeff)
            deltas.append(pur - base_pur)
            rows.append({
                "parameter": name,
                "perturbation_pct": sign * 5,
                "value": round(params[name], 4),
                "li_recovery_pct": round(rec, 4),
                "purity_pct": round(pur, 4),
                "recovery_sd_over_%d_runs" % REPEATS: round(sd, 5),
                "sensitivity_coeff": round(coeff, 5),
                "residual": model.last_residual,
            })
        print("%-12s baseline %.4f | +5%% purity %+.4f | -5%% purity %+.4f | coeff %.5f"
              % (name, BASELINE[name], deltas[0], deltas[1], statistics.fmean(coeffs)))

    with open(args.csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", args.csv, len(rows), "rows")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary = {}
    for row in rows:
        summary.setdefault(row["parameter"], []).append(row["purity_pct"] - base_pur)
    labels = sorted(summary, key=lambda k: max(abs(v) for v in summary[k]))
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for i, label in enumerate(labels):
        lo, hi = min(summary[label]), max(summary[label])
        ax.barh(i, lo, color="#c0504d")
        ax.barh(i, hi, color="#4f81bd")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("delta product purity vs baseline (percentage points), +/-5% parameter sweep")
    ax.set_title("Li2CO3 SX sensitivity tornado (%s backend, %d repeats/point)" % (args.backend, REPEATS))
    fig.tight_layout()
    fig.savefig(args.png, dpi=200)
    print("wrote", args.png)
    print("summary:", json.dumps({k: [round(v, 4) for v in vals] for k, vals in summary.items()}))


if __name__ == "__main__":
    main()
