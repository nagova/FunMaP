import os
import sys

# Ensure DLL paths are present
for p in [r"C:\Users\admin\miniforge3\Library\bin", r"C:\Users\admin\miniforge3\envs\ubermag_env\Library\bin"]:
    if os.path.exists(p) and p not in os.environ.get("PATH", ""):
        os.environ["PATH"] = p + ";" + os.environ.get("PATH", "")

import traceback


def main():
    try:
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm
        import pandas as pd

        BASE_DIR = r"c:\Users\admin\Documents\Coding Projects\FunMaP"
        SNAP_DIR = os.path.join(BASE_DIR, "results", "r1_1nm_snapshots")
        OUT_FIG_DIR = os.path.join(BASE_DIR, "results", "figures")
        ART_DIR = r"C:\Users\admin\.gemini\antigravity\brain\21f12b46-1dcb-4cc6-b36c-a2126717ebc6"

        CELL_NM = 1.0
        R_SUB = 100.0 # nm
        T0_CAP = 60.0 # nm
        THETA_MAX = 85.0 * np.pi / 180.0

        states = [
            (0.0, "Frame 1: Remanence State ($B = 0.00$ T)\nUniform radial easy-axis orientation ($M_z/M_s = +0.710$)", "frame_remanence_D200nm_1p0nm.png"),
            (-6.2, "Frame 2: Pre-Collapse State ($B = -6.20$ T)\nEquatorial rim reversed ($M_z < 0$), domain wall advanced ($M_z/M_s = +0.247$)", "frame_preswitch_D200nm_1p0nm.png"),
            (-6.3, "Frame 3: Post-Collapse State ($B = -6.30$ T)\nBarkhausen collapse of polar core ($M_z/M_s = -0.866$)", "frame_switched_D200nm_1p0nm.png"),
        ]

        for b_tgt, title, fname in states:
            npz_path = os.path.join(SNAP_DIR, f"r1_1nm_state_B{b_tgt:+.2f}T.npz")
            print(f"Loading {npz_path}...")
            data = np.load(npz_path)
            m_3d = data['m']
            norm_3d = data['norm']
            mz_val = float(data['mz'])
            
            nx, ny, nz, _ = m_3d.shape
            iy_mid = ny // 2
            cell_nm = CELL_NM
            Lx_nm = nx * cell_nm
            Lz_nm = nz * cell_nm

            x_coords = np.linspace(-Lx_nm/2 + cell_nm/2, Lx_nm/2 - cell_nm/2, nx)
            z_coords = np.linspace(cell_nm/2, Lz_nm - cell_nm/2, nz)
            X, Z = np.meshgrid(x_coords, z_coords)

            m_slice = m_3d[:, iy_mid, :, :]
            norm_slice = norm_3d[:, iy_mid, :, 0]
            mx_2d = m_slice[:, :, 0].T
            mz_2d = m_slice[:, :, 2].T
            norm_2d = norm_slice.T
            mask = norm_2d > 0.1

            # Normalize to true unit vectors
            mx_unit = np.where(mask, mx_2d / norm_2d, 0.0)
            mz_unit = np.where(mask, mz_2d / norm_2d, 0.0)


            fig = plt.figure(figsize=(11, 7.5), dpi=300)
            ax = fig.add_axes([0.08, 0.10, 0.62, 0.80])
            ax_inset = fig.add_axes([0.74, 0.50, 0.23, 0.38])
            cax = fig.add_axes([0.74, 0.15, 0.21, 0.04])

            norm_cmap = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
            mz_masked = np.ma.masked_where(~mask, mz_2d)
            im = ax.pcolormesh(X, Z, mz_masked, cmap='coolwarm', norm=norm_cmap, shading='auto')

            # Quiver with proper scale: length of arrow ~ 8 nm
            skip = 10
            X_sub = X[::skip, ::skip]
            Z_sub = Z[::skip, ::skip]
            mx_sub = mx_unit[::skip, ::skip]
            mz_sub = mz_unit[::skip, ::skip]
            mask_sub = mask[::skip, ::skip]

            ax.quiver(X_sub[mask_sub], Z_sub[mask_sub],
                      mx_sub[mask_sub], mz_sub[mask_sub],
                      color='black', scale_units='xy', angles='xy', scale=0.12,
                      width=0.0035, headwidth=3.5, headlength=4.5, headaxislength=4, pivot='mid')

            theta_sub = np.linspace(0, np.pi, 200)
            ax.plot(R_SUB * np.cos(theta_sub), R_SUB * np.sin(theta_sub), 'k--', lw=1.5, label='Substrate ($R = 100$ nm)')

            theta_cap = np.linspace(-THETA_MAX, THETA_MAX, 200)
            r_out = R_SUB + T0_CAP * np.cos(theta_cap)
            ax.plot(r_out * np.sin(theta_cap), r_out * np.cos(theta_cap), 'k-', lw=1.5, label=r'Cap surface ($t_0\cos\theta$)')

            ax.set_aspect('equal')
            ax.set_xlim(-165, 165)
            ax.set_ylim(-15, 175)
            ax.set_xlabel('Spatial Coordinate $X$ (nm)', fontsize=11, fontweight='bold')
            ax.set_ylabel('Spatial Coordinate $Z$ (nm)', fontsize=11, fontweight='bold')
            ax.set_title(title, fontsize=11.5, fontweight='bold', pad=10)
            ax.grid(True, linestyle=':', alpha=0.4)
            ax.legend(loc='lower center', fontsize=9, framealpha=0.9)

            # Inset
            r1_csv = os.path.join(BASE_DIR, "results", "loops", "R1_hemi_d200nm_cell_1p0nm.csv")
            if os.path.exists(r1_csv):
                df_l = pd.read_csv(r1_csv)
                ax_inset.plot(df_l['B_ext_T'], df_l['Mz'], 'k-', lw=1.8, label='R1 ($1.0$ nm)')
            ax_inset.plot(b_tgt, mz_val, 'ro', markersize=8, markeredgecolor='black', label=f'Current: {b_tgt:+.2f} T')
            ax_inset.axhline(0, color='gray', lw=0.6, ls=':')
            ax_inset.axvline(0, color='gray', lw=0.6, ls=':')
            ax_inset.set_xlim(-15.5, 15.5)
            ax_inset.set_ylim(-1.08, 1.08)
            ax_inset.set_xlabel(r'$\mu_0 H_z$ (T)', fontsize=9, fontweight='bold')
            ax_inset.set_ylabel(r'$M_z / M_s$', fontsize=9, fontweight='bold')
            ax_inset.set_title(r'R1 Benchmark Loop ($H_c = 6.22$ T)', fontsize=9.0, fontweight='bold')
            ax_inset.grid(True, linestyle=':', alpha=0.5)
            ax_inset.legend(loc='lower left', fontsize=7.5)

            cb = plt.colorbar(im, cax=cax, orientation='horizontal')
            cb.set_label(r'Normalized Magnetization $m_z = M_z / M_s$', fontsize=9.5, fontweight='bold')
            cb.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])

            fig.text(0.5, 0.02,
                     f"RUN R1 VERIFIED BENCHMARK (Cell = 1.0 nm): Whole Hemisphere 3D ($D=200$ nm, $t_0=60$ nm, Set A: pure $L1_0$, $0$ K)\n"
                     f"Mesh: 16.38 million cells. Inward domain wall advances from rim (85 deg); polar core switches at Hc = 6.22 T.",
                     ha='center', fontsize=9.2, fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.4', facecolor='#ebf5fb', edgecolor='#2980b9', lw=1.2))

            fig.savefig(os.path.join(OUT_FIG_DIR, fname), dpi=300)
            fig.savefig(os.path.join(ART_DIR, fname), dpi=300)
            plt.close(fig)
            print(f"Rendered {fname} successfully!")

        # Also generate a combined 3-panel publication figure!
        fig, axes = plt.subplots(1, 3, figsize=(18, 6.2), dpi=300)
        for idx, (b_tgt, title, fname) in enumerate(states):
            ax = axes[idx]
            npz_path = os.path.join(SNAP_DIR, f"r1_1nm_state_B{b_tgt:+.2f}T.npz")
            data = np.load(npz_path)
            m_3d = data['m']
            norm_3d = data['norm']
            mz_val = float(data['mz'])
            
            nx, ny, nz, _ = m_3d.shape
            iy_mid = ny // 2
            cell_nm = CELL_NM
            Lx_nm = nx * cell_nm
            Lz_nm = nz * cell_nm

            x_coords = np.linspace(-Lx_nm/2 + cell_nm/2, Lx_nm/2 - cell_nm/2, nx)
            z_coords = np.linspace(cell_nm/2, Lz_nm - cell_nm/2, nz)
            X, Z = np.meshgrid(x_coords, z_coords)

            m_slice = m_3d[:, iy_mid, :, :]
            norm_slice = norm_3d[:, iy_mid, :, 0]
            mx_2d = m_slice[:, :, 0].T
            mz_2d = m_slice[:, :, 2].T
            norm_2d = norm_slice.T
            mask = norm_2d > 0.1

            mx_unit = np.where(mask, mx_2d / norm_2d, 0.0)
            mz_unit = np.where(mask, mz_2d / norm_2d, 0.0)

            norm_cmap = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
            mz_masked = np.ma.masked_where(~mask, mz_2d)
            im = ax.pcolormesh(X, Z, mz_masked, cmap='coolwarm', norm=norm_cmap, shading='auto')

            skip = 10
            X_sub = X[::skip, ::skip]
            Z_sub = Z[::skip, ::skip]
            mx_sub = mx_unit[::skip, ::skip]
            mz_sub = mz_unit[::skip, ::skip]
            mask_sub = mask[::skip, ::skip]

            ax.quiver(X_sub[mask_sub], Z_sub[mask_sub],
                      mx_sub[mask_sub], mz_sub[mask_sub],
                      color='black', scale_units='xy', angles='xy', scale=0.12,
                      width=0.004, headwidth=3.5, headlength=4.5, headaxislength=4, pivot='mid')

            theta_sub = np.linspace(0, np.pi, 200)
            ax.plot(R_SUB * np.cos(theta_sub), R_SUB * np.sin(theta_sub), 'k--', lw=1.2)
            theta_cap = np.linspace(-THETA_MAX, THETA_MAX, 200)
            r_out = R_SUB + T0_CAP * np.cos(theta_cap)
            ax.plot(r_out * np.sin(theta_cap), r_out * np.cos(theta_cap), 'k-', lw=1.2)

            ax.set_aspect('equal')
            ax.set_xlim(-165, 165)
            ax.set_ylim(-15, 175)
            ax.set_xlabel('Spatial Coordinate $X$ (nm)', fontsize=10, fontweight='bold')
            if idx == 0:
                ax.set_ylabel('Spatial Coordinate $Z$ (nm)', fontsize=10, fontweight='bold')
            sub_title = [
                f"(a) Remanence ($B = 0.00$ T)\n$M_z/M_s = {mz_val:+.3f}$",
                f"(b) Pre-Switch ($B = -6.20$ T)\n$M_z/M_s = {mz_val:+.3f}$ (Rim inverted)",
                f"(c) Switched ($B = -6.30$ T)\n$M_z/M_s = {mz_val:+.3f}$ (Core collapsed)"
            ][idx]
            ax.set_title(sub_title, fontsize=11, fontweight='bold', pad=8)
            ax.grid(True, linestyle=':', alpha=0.4)

        fig.subplots_adjust(bottom=0.22, top=0.88, left=0.06, right=0.95, wspace=0.18)
        cbar_ax = fig.add_axes([0.25, 0.08, 0.50, 0.035])
        cb = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cb.set_label(r'Normalized Magnetization $m_z = M_z / M_s$', fontsize=10, fontweight='bold')
        cb.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])

        triptych_p1 = os.path.join(OUT_FIG_DIR, "r1_1nm_triptych_snapshots.png")
        triptych_p2 = os.path.join(ART_DIR, "r1_1nm_triptych_snapshots.png")
        fig.savefig(triptych_p1, dpi=300)
        fig.savefig(triptych_p2, dpi=300)
        plt.close(fig)
        print("Saved 3-panel triptych:", triptych_p1)
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
