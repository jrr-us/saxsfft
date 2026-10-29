# SAXS‑FFT 2D

Joaquin Ramirez-Rico, University of Seville — [mates-lab.es](https://mates-lab.es/)

**Run it online: <https://jrr-us.github.io/saxsfft/>**

A small teaching application for small‑angle scattering. It generates two‑dimensional microstructures and computes their scattering pattern as the squared modulus of the 2D Fourier transform. It shows the microstructure, the 2D pattern and the azimuthally averaged I(q).

## Running it

Use the online version above, or download this repository and open `index.html` in a recent browser (Chrome, Edge or Firefox) on Windows, Linux or macOS. Nothing needs to be installed and no internet connection is needed. You can copy the whole folder to a USB stick or put it on any web server.

The computation runs in a background Web Worker, so the page stays responsive. *Stop* aborts a long calculation.

## What it computes

The sample is an N × N grid of pixels of size *a* (nm) with periodic boundary conditions. The box size is L = N·a.

| quantity | definition |
|---|---|
| SLD map | ρ(x, y). Shapes are antialiased with 4×4 sub‑pixel sampling at their edges. |
| 2D intensity | I(**q**) = \|a² Σ Δρ(**x**) e^(−i**q**·**x**)\|² / A, with A = L² (intensity per unit area, units Δρ²·nm²) |
| q grid | q = 2πk / L, from q_min = 2π/L to q_max = π/a (Nyquist) |
| 1D profile | azimuthal average in rings of width Δq = 2π/L |
| Analytic particles | n⟨\|F(q)\|²⟩/A summed over the particles actually generated, i.e. the curve for S(q) = 1. For a disk, F = Δρ πR² · 2J₁(qR)/(qR). Ellipse and core–shell amplitudes are orientation‑averaged. |
| S(q) | I(q) divided by the analytic particle curve. Points near zeros of P(q) are hidden. |
| Porod | 2D Porod law I → 2Δρ² (P/A) q⁻³, with the interface length P measured from Σ\|∇ρ\| of the map |
| Guinier | 2D form I₀ exp(−q²Rg²/2). For a disk Rg² = R²/2. DLA uses the cluster Rg. |
| Invariant | Σ I (Δq/2π)² = ⟨Δρ²⟩ exactly (Parseval). The info panel also shows how much of it the 1D integral up to q_max recovers. |

**Resolution smearing** (on by default) mimics a real instrument, where each detector pixel collects a small range of q (beam size, divergence, pixel size). The 2D intensity is convolved with a Gaussian in q of the chosen FWHM, given in units of Δq = 2π/L (default 3 Δq), before the radial average. The analytic curve is smeared with the same 2D Gaussian (exact radial kernel), and only the smeared version is shown. Smearing fills in the form‑factor minima, averages out speckle in the 2D pattern and leaves the invariant unchanged. Seen from real space, it multiplies the density autocorrelation γ(r) by exp(−r²σ_q²/2): the instrument cannot see correlations much longer than 1/σ_q.

The **“Include pixel blur”** option multiplies the analytic curve by the pixel transfer function ⟨sinc²(qₓa/2) sinc²(q_ya/2)⟩. Rasterising with area coverage convolves the structure with one pixel. Without the correction, the FFT data fall below the analytic curve near q_max.

**Display.** The 2D pattern can be shown as log₁₀ I (default; *2D range* sets the number of decades) or as linear I. In both, *2D top* lowers the top of the colour scale (by decades in log, as a fraction of the maximum in linear) to bring out weaker features. A linear map shows how strongly the intensity is concentrated at low q, which the log map hides.

## Structures

- **Single disk**: textbook form factor. Minima at qR = 3.83, 7.02, 10.17.
- **Dilute disks**: random positions without overlap. Optional polydispersity (Gaussian or log‑normal) smears out the minima. *Overlaps allowed* gives a Poisson (ideal‑gas) arrangement.
- **Hard‑disk fluid (dense)**: φ up to 0.72. Overlaps are removed at full size by relaxation, then Metropolis hard‑disk Monte Carlo equilibrates the fluid (the *MC sweeps* setting). Shows the correlation peak and S(0) < 1.
- **Ellipses**: random or aligned orientation (angle and spread) give isotropic or anisotropic patterns. Exact non‑overlap uses the Perram–Wertheim contact function.
- **Core–shell disks**: three SLD levels. Try contrast matching: shell = matrix, or core = shell.
- **Two‑phase random (spinodal)**: a thresholded Gaussian random field with a ring spectrum (broad peak at 2π/d₀) or a Gaussian spectrum.
- **Fractal aggregates (DLA)**: off‑lattice diffusion‑limited aggregation. The info panel reports the mass–radius dimension. With 10³ particles per cluster the power‑law range is short, so the apparent slope is smaller than D_f ≈ 1.71: a nice discussion point.

## Things to try in class

1. Single disk: read R from the first minimum (qR = 3.83) and from the Guinier slope. Switch on the Porod line.
2. Dilute → dense hard disks at the same R: watch I(q→0) drop and the peak appear near 2π/(2R). Open the S(q) tab.
3. Polydispersity 0 → 0.2: the form‑factor minima vanish.
4. Tick *Invert phases*. I(q) is unchanged for q > 0 (Babinet's principle).
5. Untick *Subtract mean ρ* and enable zero padding: the q = 0 peak leaks into neighbouring q (truncation and sinc ripples). With smearing on, it spreads like a direct beam: the reason for a beamstop.
6. Single disk: compare with *Resolution smearing* off and on, then increase the FWHM step by step. The minima fill in, while the Guinier and Porod regions hardly change.
7. Aligned ellipses: the pattern is elongated perpendicular to the long axis (reciprocal space).
8. Make the box smaller (N = 256): q_min rises and the Guinier region is lost. Make the pixel bigger: q_max falls.

## Export

- PNG of each panel.
- CSV of q, I, the analytic particle curve, S(q) and the pixel count per ring. The parameters are written in the header.
- *Save params* / *Load params* store a JSON file, so a figure can be reproduced exactly. All randomness is seeded.

## Files

| file | content |
|---|---|
| `index.html`, `css/style.css` | page layout |
| `js/fft.js` | radix‑2 FFT, 2D FFT, fftshift |
| `js/rng.js` | seeded random numbers, size distributions |
| `js/structures.js` | structure generators, packing, DLA |
| `js/raster.js` | antialiased rasterisation into the SLD map |
| `js/analysis.js` | intensity, radial average, analytic curves, invariant, Porod |
| `js/worker.js` | Web Worker built from the files above (works from `file://`) |
| `js/plot.js`, `js/colormaps.js` | canvas plotting |
| `js/main.js` | user interface |

The computational files are plain JavaScript with no dependencies, so they can also be loaded in Node.js for scripting.

## Licence

MIT, see [LICENSE](LICENSE). © 2026 Joaquin Ramirez-Rico.
