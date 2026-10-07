---
marp: true
theme: default
paginate: true
header: '**FunMaP v2** | Evolution of Micromagnetic Modeling in FePt Janus Caps'
footer: 'Group Meeting Presentation — Physical Rationale & Resolution of Brown’s Paradox'
style: |
  section {
    font-family: 'Segoe UI', -apple-system, BlinkMacSystemFont, Roboto, sans-serif;
    padding: 32px 48px;
    background-color: #0f172a;
    color: #f8fafc;
  }
  h1 {
    color: #ffffff;
    font-size: 24px;
    margin-bottom: 4px;
    border-bottom: 2px solid #14b8a6;
    padding-bottom: 6px;
  }
  h2 {
    color: #14b8a6;
    font-size: 14px;
    text-transform: uppercase;
    letter-spacing: 1px;
    margin-bottom: 8px;
  }
  .subtitle {
    color: #94a3b8;
    font-size: 13px;
    margin-bottom: 16px;
  }
  .columns-3 {
    display: grid;
    grid-template-columns: 1fr 1fr 1fr;
    gap: 16px;
  }
  .columns-2 {
    display: grid;
    grid-template-columns: 1.15fr 0.85fr;
    gap: 20px;
  }
  .card {
    background: #1e293b;
    border: 1px solid #475569;
    border-radius: 8px;
    padding: 12px 14px;
    font-size: 11px;
    line-height: 1.45;
  }
  .card-red { border-top: 3px solid #ef4444; }
  .card-blue { border-top: 3px solid #3b82f6; }
  .card-orange { border-top: 3px solid #f97316; }
  .card-teal { border-top: 3px solid #14b8a6; }
  .card-green { border-top: 3px solid #22c55e; }
  .card-title {
    font-weight: 700;
    font-size: 13px;
    margin-bottom: 8px;
    color: #ffffff;
  }
  .badge {
    display: inline-block;
    padding: 2px 6px;
    border-radius: 4px;
    font-size: 9px;
    font-weight: 700;
    text-transform: uppercase;
    float: right;
  }
  .badge-red { background: rgba(239, 68, 68, 0.25); color: #fca5a5; }
  .badge-blue { background: rgba(59, 130, 246, 0.25); color: #93c5fd; }
  .badge-orange { background: rgba(249, 115, 22, 0.25); color: #fdba74; }
  .badge-teal { background: rgba(20, 184, 166, 0.25); color: #5eead4; }
  .badge-green { background: rgba(34, 197, 94, 0.25); color: #86efac; }
  .highlight-box {
    background: #0f172a;
    border-left: 3px solid #14b8a6;
    padding: 6px 10px;
    margin-top: 8px;
    font-weight: 600;
    color: #38bdf8;
    font-size: 11px;
  }
  .banner {
    background: #1e293b;
    border: 1px solid #14b8a6;
    border-radius: 8px;
    padding: 10px 14px;
    margin-top: 14px;
    font-size: 11.5px;
    color: #f1f5f9;
  }
  table {
    width: 100%;
    border-collapse: collapse;
    font-size: 10.5px;
    margin-top: 6px;
  }
  th {
    background: #0f172a;
    color: #94a3b8;
    text-align: left;
    padding: 4px 6px;
    border-bottom: 1px solid #475569;
  }
  td {
    padding: 4px 6px;
    border-bottom: 1px solid rgba(255,255,255,0.06);
    color: #e2e8f0;
  }
  ul {
    margin: 0;
    padding-left: 14px;
  }
  li {
    margin-bottom: 4px;
  }
---

## Diagnostic Campaign & Brown's Paradox
# 1. Why Legacy Models Failed & What Diagnostic Controls Proved
<div class="subtitle">Tracing the transition from coarse-grid discretization artifacts to rigorous single-crystal & disorder bounds</div>

<div class="columns-3">

<div class="card card-red">
<span class="badge badge-red">Identified Artifacts</span>
<div class="card-title">Legacy v1 Model Breakdown</div>

- **Discretization Breakdown:** Grid cells $\Delta x = 18\text{--}92\text{ nm} \gg \ell_{\mathrm{ex}} \approx 2.1\text{ nm}$. Underestimated inter-cell exchange coupling, causing artificial spin decoupling and unphysical remanence collapse ($M_r/M_s \to 0.40$).
- **Geometric Error:** Assumed uniform $t = 60\text{ nm}$ shell with blunt $90^\circ$ rims, placing unphysical magnetic volume at the equator.
- **Bisection Solver Trap:** Tested negative branch without proper sweep history, producing an artificial plateau at $\mu_0 H_c = 4.71\text{ T}$ and a $12.2\text{ T}$ numerical ceiling.

<div class="highlight-box" style="border-color: #ef4444; color: #fca5a5;">
Key Artifact: 4.71 T plateau & 12.2 T ceiling were purely numerical.
</div>
</div>

<div class="card card-blue">
<span class="badge badge-blue">Pure Geometry Control</span>
<div class="card-title">Continuum Benchmark (R1)</div>

- **Physical Geometry:** Sub-exchange mesh ($\Delta x = 1.0\text{ nm} < \ell_{\mathrm{ex}}$), ballistic taper $t(\theta) = t_0\cos\theta$, radial easy axes, 100% pure $L1_0$ phase ($0\%$ A1).
- **Monotonic Sweeps:** Eliminated solver trapping via descending field sweep from $+7\text{ T}$.
- **Brown's Paradox Confirmed:** Single-crystal cap yields $\mu_0 H_c = \mathbf{6.22\text{ T}}$ ($M_r/M_s = 0.71$) via rim nucleation. Curvature alone lowers $H_c$ from Stoner-Wohlfarth ceiling ($13.2\text{ T}$), but **cannot reach experimental 1.13 T**.

<div class="highlight-box" style="border-color: #3b82f6; color: #93c5fd;">
Finding: Curvature alone ≠ 1.13 T (Discrepancy: ΔHc = +5.09 T).
</div>
</div>

<div class="card card-orange">
<span class="badge badge-orange">The Remanence Dilemma</span>
<div class="card-title">Uniform Disorder Control (M1)</div>

- **Homogeneous Soft Phase:** Evaluated caps with uniform A1 fraction ($f_{A1} = 0\text{--}66\%$) at $\Delta x = 1.3\text{ nm}$. Soft grains lower coercivity at $-7.1\text{ T}/f_{A1}$.
- **The Remanence Dilemma:** Lowering coercivity to $1.13\text{ T}$ requires $f_{A1} > 75\%$, which collapses remanence below $M_r/M_s < 0.30$.
- **Reviewer Critique Validated:** Homogeneous disorder is fundamentally incapable of simultaneously satisfying $\mu_0 H_c \approx 1.13\text{ T}$ AND high remanence $M_r/M_s \approx 0.64$.

<div class="highlight-box" style="border-color: #f97316; color: #fdba74;">
Takeaway: Uniform disorder cannot fit both Hc and Mr simultaneously.
</div>
</div>

</div>

<div class="banner">
⚖️ <strong>Fundamental Physics Conclusion:</strong> Geometry alone overestimates coercivity ($6.22\text{ T} \gg 1.13\text{ T}$), while uniform chemical disorder destroys remanence ($M_r/M_s < 0.30$). Resolving Brown's paradox requires an <strong>inhomogeneous order gradient</strong> naturally emerging from thin-film deposition kinetics.
</div>

---

## Physical Solution & Experimental Reality
# 2. Multiscale Order-Graded Architecture & Resolution of Brown’s Paradox
<div class="subtitle">Coupling line-of-sight deposition, thickness-dependent kinetics, and thermal activation</div>

<div class="columns-2">

<div class="card card-teal">
<span class="badge badge-teal">Multiscale Model (M2/M4)</span>
<div class="card-title">The Four Physical Pillars of Model v2</div>

- **Pillar I (Line-of-Sight Flux):** Cosine ballistic taper $t(\theta) = t_0\cos\theta$ ($60\text{ nm}$ at pole $\to 0\text{ nm}$ at equator).
- **Pillar II (Ordering Kinetics Gradient):** Interdiffusion and ordering scale with thickness: $S(\theta) = 0.90\sqrt{\cos\theta}$. Dewetted rim ($t < 15\text{ nm}$) acts as a soft nucleation pad ($S \to 0$), initiating reversal at $-0.3\text{ T}$.
- **Pillar III (Sub-Exchange Discretization):** Mesh cell $\Delta x = 1.3\text{ nm} < \ell_{\mathrm{ex}} \approx 2.1\text{ nm}$, strictly resolving exchange spring and domain wall boundaries ($\delta_w \approx 3.5\text{ nm}$).
- **Pillar IV (300 K Sharrock Thermal Scaling):** Thermal activation over SQUID measurement window ($\tau \approx 10\text{ s}$) reduces coercive field by $19.4\%$ for nanoscale grains ($d_g \approx 15\text{ nm}$).

<div style="background: rgba(15,23,42,0.6); padding: 8px 10px; border-radius: 6px; margin-top: 8px;">
<strong>Two-Zone Cooperative Mechanism:</strong><br>
• <strong>Polar Crown ($\theta \le 60^\circ$):</strong> $75\%$ magnetic volume, $t > 30\text{ nm}$, $S \ge 0.64 \implies$ locks remanence ($M_r/M_s = 0.637$).<br>
• <strong>Equatorial Rim ($\theta > 70^\circ$):</strong> $t < 15\text{ nm}$, disordered A1 $\implies$ triggers early nucleation without collapsing crown remanence.
</div>
</div>

<div class="card card-green">
<span class="badge badge-green">Exact Fit</span>
<div class="card-title">Quantitative Match to Experiment</div>

<table>
<thead>
<tr><th>Parameter</th><th>Model v2 (M2/M4)</th><th>Exp. SQUID (3–10 µm)</th></tr>
</thead>
<tbody>
<tr><td><strong>Coercivity ($\mu_0 H_c$)</strong></td><td style="color:#22c55e; font-weight:bold;">1.127 T (at 300 K)</td><td><strong>1.13 ± 0.04 T</strong> (Pooled mean)</td></tr>
<tr><td><strong>Remanence ($M_r/M_s$)</strong></td><td style="color:#22c55e; font-weight:bold;">0.637</td><td><strong>0.58 – 0.65</strong></td></tr>
<tr><td><strong>Diameter Scaling</strong></td><td>Invariant ($\Delta M_r < 0.01$)</td><td>Null ($p = 0.119$, &lt;13% var.)</td></tr>
<tr><td><strong>Zero-Field Stability</strong></td><td>Locked single-domain</td><td>Single vortex / uniform cap</td></tr>
</tbody>
</table>

<div class="card-title" style="margin-top: 12px; font-size: 12px; color: #94a3b8;">Physical Origins of Diameter Invariance:</div>

- **Scale Decoupling:** Exchange length $\ell_{\mathrm{ex}} \approx 2\text{ nm} \ll R = 1.5\text{--}5\ \mu\text{m}$ ($\ell_{\mathrm{ex}}/R \sim 10^{-3}$).
- **Exchange Spring Pinning:** Nucleation initiated at the rim is governed by the local thickness gradient, not sphere diameter.
- **Interparticle Necks:** Physical bridge percolation ($\eta \approx 0.90$) stabilizes adjacent caps against collective reversal.

<div class="highlight-box" style="border-color: #22c55e; color: #86efac;">
Result: Comprehensive, parameter-consistent explanation of experimental FePt cap physics.
</div>
</div>

</div>
