# rgd-kernel

**Relational Gravitational Dynamics — Validation Pipeline**

Public code for testing the RGD effective gravitational kernel against four
independent observational regimes, as presented in:

> Maia, K. (2026). *A Cosmological-Background Effective Gravitational Kernel:
> Unifying disk galaxies, galaxy groups, galaxy clusters, and gravitational
> lensing under a single effective law.* Relational Gravitation Program, v18.

---

## The kernel

$$K_{\rm eff}(g, z;\, w, \xi) = 1 + \frac{\lambda\, e^{-g/g^\dagger(z)}}{\left[1 - e^{-g/g^\dagger(z)}\right]^{n(w,\xi)}}$$

| Symbol | Definition | Value / source |
|---|---|---|
| $g^\dagger(z)$ | $c\,H(z)/(2\pi)$ | cosmological scale, no free parameters |
| $\lambda$ | amplitude | 0.999 (fitted on SPARC) |
| $n(w,\xi)$ | $(1/2 + A\,\xi)\,w^{-\alpha(\xi)}$ | transition exponent |
| $A$ | exponent-law coefficient | 0.4068 ± 0.005 (59 systems) |
| $\xi$ | $g_{500}/g^\dagger(z)$ | depth in cosmological regime |
| $w$ | $r/R_{500}$ | normalized radius |
| $\alpha(\xi)$ | $0.355\,(1-e^{-\xi/0.052})$ | radial modulation |

**Disk-galaxy limit** ($\xi \to 0$): $n = 1/2$, derived — not fitted — from the
flat-rotation-curve condition. The Baryonic Tully-Fisher Relation follows
analytically as $v_f^4 = \lambda^2 g^\dagger G M_b$.

---

## Performance summary

| Sample | N | Observable | RMS |
|---|---|---|---|
| SPARC disk galaxies | 171 | Rotation curves | 0.299 dex |
| Sun+2009 groups | 23 | Chandra HSE | 0.034 dex |
| X-COP clusters (w=1) | 13 | X-ray + SZ | 0.015 dex |
| CLASH lensing | 23 | Weak lensing | 0.005 dex (in-sample) |
| X-COP radial profiles | 13 × 3 | X-ray at 3 radii | 0.040 dex |

LOO-CLASH predictive RMS: **0.006 dex**.  
Radial test vs MOND: **6.3× better** (7.8× in the fully predictive subset).

---

## Repository contents

```
rgd_kernel.py            Main pipeline — kernel, tests, calibration, plots
sparc_gb.npy             g_bar for SPARC (175 galaxies, 3391 points) [m/s²]
sparc_go.npy             g_obs for SPARC [m/s²]
empirical_xi.npy         ξ values for groups + X-COP (59 systems)
empirical_n.npy          Empirical n values for groups + X-COP
empirical_labels.npy     Labels: 'G' = group, 'C' = X-COP cluster
sun2009_processed.txt    Sun+2009 group table (Chandra HSE)
```

X-COP (Ettori+2019) and CLASH (Umetsu+2016) data are embedded directly in
`rgd_kernel.py` — no additional files required for those samples.

---

## Installation

```bash
git clone https://github.com/Maia-KarlosMarden/rgd-kernel.git
cd rgd-kernel
pip install numpy scipy matplotlib pandas   # matplotlib and pandas optional
```

No package installation step — the pipeline runs as a single script.

---

## Usage

```bash
# Run all five tests and print the summary table
python rgd_kernel.py

# Run a single sample
python rgd_kernel.py --sample sparc
python rgd_kernel.py --sample groups
python rgd_kernel.py --sample xcop
python rgd_kernel.py --sample xcop_radial
python rgd_kernel.py --sample clash

# Generate the six-panel validation figure
python rgd_kernel.py --plot

# Export all residuals to CSV
python rgd_kernel.py --export results.csv

# Re-calibrate A from the embedded data (reproduces A = 0.4068)
python rgd_kernel.py --calibrate

# Calibrate A without using CLASH (tests out-of-sample prediction)
python rgd_kernel.py --calibrate --no-clash-cal
```

### Use as a library

```python
from rgd_kernel import kernel, g_dagger, n_exponent, g_pred_rgd, rms_dex
import numpy as np

# Kernel at a single point
g_bar = 1e-11          # m/s²  (deep-MOND regime)
z     = 0.5
xi    = 0.3            # typical for a galaxy group
n     = n_exponent(xi, w=1.0)
K     = kernel(g_bar, z, n)
print(f"n = {n:.3f},  K_eff = {K:.2f}")

# SPARC test from pre-computed arrays
from rgd_kernel import test_sparc
r = test_sparc(datadir=".")
print(f"SPARC RMS = {r['rgd']['rms']:.4f} dex")

# Exact 171-galaxy reproduction (requires Rotmod_LTG/ from SPARC website)
r = test_sparc(rotmod_dir="/path/to/Rotmod_LTG/")
```

---

## Reproducing the paper results exactly

### SPARC (171 galaxies)

The pre-computed arrays `sparc_gb.npy` / `sparc_go.npy` include all 175
galaxies. Four are excluded in the paper on physical grounds (see §3.2):

```python
SPARC_EXCLUDED = ['CamB', 'NGC4389', 'UGC02455', 'UGC11557']
```

The point-level filter `g_obs > 0.5 × g_bar` removes most of their anomalous
points; the RMS difference from including all 175 is < 0.001 dex. For exact
reproduction, download `Rotmod_LTG/` from
[astroweb.cwru.edu/SPARC](http://astroweb.cwru.edu/SPARC/) and pass
`--rotmod_dir /path/to/Rotmod_LTG/`.

### X-COP and CLASH

Data embedded in `rgd_kernel.py` (`xcop_data()`, `xcop_radial_data()`,
`clash_data()`). No external download required. Original sources:

- Ettori et al. (2019), A&A 621 A39 — X-COP hydrostatic masses
- Ghirardini et al. (2019), A&A 621 A41 — X-COP gas fractions
- Umetsu et al. (2016), ApJ 821 116 — CLASH lensing masses

### Sun+2009 groups

Loaded from `sun2009_processed.txt`. Original source:
Sun et al. (2009), ApJ 693 1142.

---

## Key physical notes

**g†(z) evolution.** The cosmological scale $g^\dagger(z) = c H(z)/(2\pi)$
grows with redshift: +27% at z=0.40 (mean CLASH), +35% at z=0.55, +65% at
z=0.89. This evolution is empirically detected in CLASH at 4.1σ (permutation
test) and drives the residual correlation r=0.999.

**α(ξ) calibration status.** The radial modulation $\alpha(\xi)$ was calibrated
on 8 X-COP clusters in an earlier pipeline version (2 free parameters). Of the
39 points in the radial test, 15 (5 clusters) are genuinely predictive;
RMS = 0.031 dex and 7.8× improvement over MOND for that subset. Full
recalibration against eROSITA DR1 profiles is planned for Phase II.

**Equilibrium domain.** The kernel is valid for systems in circular or
hydrostatic equilibrium. Active mergers (e.g. El Gordo, z=0.87) lie outside
this domain; their lensing mass includes kinetic contributions not modelled
by the kernel.

---

## Citation

```bibtex
@misc{maia2026rgd,
  author  = {Maia, Karlos},
  title   = {A Cosmological-Background Effective Gravitational Kernel},
  year    = {2026},
  note    = {Relational Gravitation Program},
  url     = {https://github.com/Maia-KarlosMarden/rgd-kernel}
}
```

---

## License

MIT. Attribution required — please cite the paper above if you use this code
in published work.

---

## Update: v5 (October 2026)

A second paper extends this pipeline: the radial modulation `w^(-α(ξ))` is replaced by a one-parameter scale shift
`s = [Δ(r)/Δ500]^p` (X-COP 39-point rms 0.037 → 0.016–0.021 dex; works at r2500 where the v4 law diverges), g†(z) is tested
on Allen and Mantz f_gas, and groups, lensing masses, the high-z BTFR and a cosmological damper are added. Several
tensions are documented rather than hidden (f_bar–mass trend at r2500, ~20% normalisation offset between samples,
mild preference of z ≈ 1–2 galaxies for a weaker g†(z) evolution). See [`v5/README.md`](v5/README.md).
The v4 kernel table above is unchanged and the original pipeline remains in `rgd_pipeline.zip`.
