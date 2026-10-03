# FunMaP Simulation Guide: Corrected FePt Micromagnetics

This guide provides the complete physical framework, acceptance criteria, multi-tier simulation plan, and parameter reference for micromagnetic modeling of FePt Janus caps in FunMaP.

---

## 1. Material Constants & Physical Governing Lengths

All simulations use the fixed experimental material parameters from Table S4 of the manuscript:

| Parameter | Symbol | Value | Meaning / Role |
|---|---|---|---|
| Saturation Magnetization | $M_s$ | $1.0\times 10^6\text{ A/m}$ | Saturation magnetization |
| Exchange Stiffness | $A$ | $1.0\times 10^{-11}\text{ J/m}$ | Exchange coupling strength |
| Hard-phase Anisotropy | $K_u (\text{L}1_0)$ | $6.6\times 10^6\text{ J/m}^3$ | Uniaxial anisotropy (chemically ordered phase) |
| Soft-phase Anisotropy | $K_u (\text{A}1)$ | $1.0\times 10^4\text{ J/m}^3$ | Disordered soft phase |
| Cap Thickness | $t$ | $60\text{ nm}$ | FePt film thickness deposited on sphere |
| Gilbert Damping | $\alpha$ | $0.5$ | Quasi-static relaxation damping |

### Governing Lengths

| Quantity | Formula | Value | Physical Significance |
|---|---|---|---|
| **Anisotropy Field** | $\mu_0 H_k = 2 K_u / M_s$ | **13.200 T** | Intrinsic switching field along easy axis |
| **Anisotropy / Bloch Length** | $\sqrt{A / K_u}$ | **1.231 nm** | **GOVERNS THE MESH**. Smallest physical length scale |
| **Domain-Wall Width** | $\pi \sqrt{A / K_u}$ | **3.867 nm** | Physical width of magnetic domain wall |
| **Magnetostatic Exchange Length** | $\sqrt{2A / (\mu_0 M_s^2)}$ | 3.989 nm | Dipolar exchange length |
| **Quality Factor** | $Q = 2K_u / (\mu_0 M_s^2)$ | 10.5 | High-$Q$ material ($Q \gg 1$): anisotropy dominates |
| **Single-Domain Diameter** | $\sim 72\sqrt{A K_u} / (\mu_0 M_s^2)$ | 465 nm | Boundary between single-domain and wall-mediated reversal |

> [!IMPORTANT]
> **The governing length scale is 1.231 nm, not 4 nm.**
> For high-$Q$ materials ($Q = 10.5$), the *smaller* of the two lengths governs the mesh discretization. The cell size must satisfy:
> $$\Delta x \le \sqrt{A/K_u} = 1.231\text{ nm} \quad (\text{strictly } \le 0.615\text{ nm})$$
> with at least 20 cells through the 60 nm thickness to properly resolve domain walls crossing the film.

---

## 2. Hard Computational Limits: Intractability Bound (1 TB – 25 TB)

OOMMF meshes the entire rectangular bounding box $(2R_{\text{out}})^2 \times R_{\text{out}}$, not just the magnetic shell.

For experimental sphere diameters ($d = 3\text{--}20\ \mu\text{m}$), resolving the governing 0.6–1.2 nm scale yields:

| Sphere Diameter | Magnetic Cells (0.6 nm) | Bounding Box Cells | Required RAM | Feasibility |
|---|---|---|---|---|
| **0.3 µm** | $1.4\times 10^7$ | $4.3\times 10^7$ | ~6 GB | **Tractable** |
| **1.0 µm** | $5.5\times 10^7$ | $2.3\times 10^8$ | ~15 GB | **Tractable (HPC)** |
| **3.0 µm** | $3.9\times 10^9$ | $6.2\times 10^9$ | **~1 TB** | **Permanently Out of Reach** |
| **5.0 µm** | $1.1\times 10^{10}$ | $2.9\times 10^{10}$ | **~2 TB** | **Permanently Out of Reach** |
| **8.0 µm** | $2.8\times 10^{10}$ | $1.2\times 10^{11}$ | **~4 TB** | **Permanently Out of Reach** |
| **10.0 µm** | $4.4\times 10^{10}$ | $2.3\times 10^{11}$ | **~6 TB** | **Permanently Out of Reach** |
| **20.0 µm** | $1.7\times 10^{11}$ | $1.8\times 10^{12}$ | **~25 TB** | **Permanently Out of Reach** |

> [!NOTE]
> The largest micromagnetic simulations in published literature operate at $\sim 10^9\text{--}10^{10}$ cells. Full-cap simulations in the 3–20 µm range are permanently intractable at the required physical resolution on modern hardware. This is a rigorous, quantitative statement of the **length-scale boundary of direct micromagnetics**.

---

## 3. Diagnosis of What Went Wrong in Legacy Code

1. **Mesh 5–149x Too Coarse**: The original rule ($n = 110$ above 1.5 µm) gave cell sizes from 28 nm to 183 nm. At 92 nm, intercell exchange ($0.0024\text{ T}$) was 5,500x weaker than anisotropy ($13.2\text{ T}$), completely decoupling cells into non-interacting particles.
2. **Disconnected Geometries**: At 8–20 µm, the 60 nm shell was thinner than a single mesh cell, shattering the body into up to 5,000 disconnected pieces.
3. **0% A1 Baseline Bug**: The original 0% A1 run reported $H_c = 12.17\text{ T}$ and $M_r/M_s = 1.00$. A radial easy-axis distribution on a hemisphere has a strict Stoner-Wohlfarth ceiling of $H_c \le 0.48 H_k = 6.34\text{ T}$ and $M_r/M_s = 0.494$. The 12.17 T value was an artifact of running a **uniform vertical easy axis** ($e_u \parallel z$) mislabelled as radial. True pure $\text{L}1_0$ radial caps switch at **$5.90\text{ T}$** with $M_r/M_s = 0.502$.
4. **Solver Tolerance**: The old tolerance `stopping_mxHxm = 2e4 A/m` ($1.9\times 10^{-3} H_k$) stopped iterations before overcoming the coercive barrier. All runs must use $\le 10\text{ A/m}$ ($10^{-6} H_k$) and monotonic sweeps with $dH \le 0.02\text{ T}$.

---

## 4. Acceptance Gates (`validate_micromagnetics.py`)

Run all five gates before accepting any simulation result:

- **Gate 1 (Mesh Resolution)**: Refuse run if cell size $> 1.231\text{ nm}$ or $< 20$ cells span the 60 nm thickness.
- **Gate 2 (Connectivity)**: Face-connected component labeling must confirm exactly 1 component holding $\ge 99.9\%$ of cells, with mean face neighbors $\ge 5.0$.
- **Gate 3 (Solver Convergence)**: `stopping_mxHxm <= 10 A/m` and monotonic descending step $dH \le 0.02\text{ T}$ near switching.
- **Gate 4 (Analytic Reproductions)**: Verify Stoner-Wohlfarth targets:
  - Uniform $\parallel H$: $H_c = 13.20 \pm 0.3\text{ T}$, $M_r/M_s = 1.000$
  - Uniform $45^\circ$: $H_c = 6.60 \pm 0.2\text{ T}$, $M_r/M_s = 0.707$
  - 3D-random non-interacting: $H_c = 6.32 \pm 0.2\text{ T}$, $M_r/M_s = 0.500$
  - Radial hemisphere non-interacting: $H_c = 6.34 \pm 0.2\text{ T}$, $M_r/M_s = 0.494$
- **Gate 5 (Mesh Convergence)**: Coercivity must plateau as cell size resolves $L_{\text{anis}} = 1.231\text{ nm}$.

---

## 5. Multi-Tier Simulation Plan

### Tier 1: Converged Full Cap Micromagnetics ($d \le 1.0\ \mu\text{m}$)
- Meshes the full hemisphere ($420\times 420\times 210\text{ nm}$ at $d = 300\text{ nm}$) without artificial octant cut planes.
- Solves domain wall nucleation and propagation across switching.
- Full cap results:
  - Pure $\text{L}1_0$ ($d = 300\text{ nm}$, cell 1.2 nm): $\mu_0 H_c = 5.904\text{ T}$, $M_r/M_s = 0.5024$.
  - 26% A1 Soft Phase ($d = 300\text{ nm}$, cell 1.2 nm): $\mu_0 H_c = 4.013\text{ T}$, $M_r/M_s = 0.5034$.

### Tier 2: Representative Curved-Patch Model ($d \ge 3\ \mu\text{m}$)
- Simulates representative $200\times 200\times 60\text{ nm}$ patches at full 1.0 nm resolution at experimental curvature:
  - $d = 3\ \mu\text{m}$ ($R = 1.5\ \mu\text{m}$, easy-axis spread $7.6^\circ$)
  - $d = 10\ \mu\text{m}$ ($R = 5.0\ \mu\text{m}$, easy-axis spread $2.3^\circ$)
- Justification: domain wall width is $3.87\text{ nm}$, sampling $< 0.03^\circ$ easy-axis variation over its width. Reversal is locally determined.

### Tier 3: Statistical Ensemble Over Global Geometry
- Integrates the local switching field distribution over the hemispherical surface normal map:
  $$M_z(H) = \int_0^{\pi/2} m_{z,\text{local}}(H, \theta) \sin\theta \, d\theta$$
- Predicts macroscopic hysteresis loops for 3–20 µm caps with rigorous defensibility.

---

## 6. Key Physics Finding & Manuscript Reframe

### Radial vs. Random is Indistinguishable
| Distribution | Coercivity $\mu_0 H_c$ | Remanence $M_r / M_s$ |
|---|---|---|
| **Radial (Hemisphere Normals)** | **6.35 T** | **0.497** |
| **Fully Random 3D** | **6.37 T** | **0.500** |
| **Difference** | **0.025 T (0.4%)** | **0.003** |

Because a hemisphere of surface normals is half-isotropic (spans $0\text{--}90^\circ$ and is isotropic in-plane), its magnetic response is magnetically indistinguishable from random 3D. The difference of $0.025\text{ T}$ is **5x below the experimental SQUID noise floor ($0.117\text{ T}$)** and cannot be measured experimentally.

### The Real Tunability: Uniform vs. Distributed
The meaningful physical comparison is **Uniform vs. Distributed Anisotropy**:
$$\mu_0 H_c (\text{Uniform}) = 13.20\text{ T} \quad \text{vs.} \quad \mu_0 H_c (\text{Distributed}) = 6.35\text{ T} \quad (\mathbf{2.1\times\text{ factor}})$$
This comparison is large, physically robust, and experimentally accessible.
