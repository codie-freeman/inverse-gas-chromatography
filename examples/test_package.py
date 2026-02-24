#!/usr/bin/env python3
"""Test IGC-SEA Package - All modules in one script."""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Import our package
from igcsea.parsing import parse_igc_csv
from igcsea.analysis.dorris_gray import prepare_alkane_data, fit_dorris_gray
from igcsea.analysis.acid_base import calculate_acid_base_params
from igcsea.analysis.surface_energy import calculate_surface_energy_profile
from igcsea.analysis.manual_dispersive import (
    calculate_dispersive_from_injections,
    prepare_alkane_data_from_injections,
    validate_against_sms,
)
from igcsea.core.constants import NA, A_CH2, R

print("=" * 60)
print("Testing IGC-SEA Package")
print("=" * 60)

# 1. Parse CSV
print("\n1. PARSING CSV")
print("-" * 60)
result = parse_igc_csv("data/examples/sample_igc_export.csv")
print(f"Parsed from: {result.source_path.name}")
print(f"  - Injection items: {result.injection_items.shape}")
if result.free_energy is not None:
    print(f"  - Free energy table: {result.free_energy.shape}")
else:
    print("  - Free energy table: not present")
if result.dispersive_surface_energy is not None:
    print(f"  - Dispersive SE table: {result.dispersive_surface_energy.shape}")
else:
    print("  - Dispersive SE table: not present")

print(f"\n  Key constants in use:")
print(f"  R  = {R} J/(mol·K)")
print(f"  NA = {NA:.8e} mol⁻¹")
print(f"  a_CH2 = {A_CH2:.15e} m²")

# 2. Dorris-Gray from SMS Free Energy table (existing pipeline)
print("\n2. DORRIS-GRAY (SMS FREE ENERGY TABLE)")
print("-" * 60)

if result.free_energy is not None:
    alkanes = prepare_alkane_data(result)
    print(f"Alkanes found: {sorted(alkanes['Solvent'].unique())}")
    fits = fit_dorris_gray(alkanes)
    print(f"\n  Linear fits (RTlnV vs Carbon Number):")
    print(fits.to_string(index=False))
else:
    print("  Skipped — free_energy table not present in this CSV.")
    alkanes = None
    fits = None

# 3. Manual Dorris-Gray from raw injection data (new independent pipeline)
print("\n3. MANUAL DORRIS-GRAY (FROM RAW INJECTION DATA)")
print("-" * 60)

_alkane_range = ["HEPTANE", "OCTANE", "NONANE", "DECANE"]
gd_com = calculate_dispersive_from_injections(result, retention_type="COM", alkanes=_alkane_range)
gd_max = calculate_dispersive_from_injections(result, retention_type="MAX", alkanes=_alkane_range)

print("  COM results:")
print(gd_com.to_string(index=False))
print("\n  MAX results:")
print(gd_max.to_string(index=False))

# 4. RTlnV vs Carbon Number — manual interpolated vs SMS
print("\n4. RTlnV vs CARBON NUMBER (MANUAL INTERPOLATED)")
print("-" * 60)

output_dir = Path(__file__).parent.parent / "outputs"
output_dir.mkdir(exist_ok=True)

alkane_inj_com = prepare_alkane_data_from_injections(result, retention_type="COM", alkanes=_alkane_range)
fits_inj_com = fit_dorris_gray(alkane_inj_com)

fig_rtlnv, axes_rtlnv = plt.subplots(2, 3, figsize=(12, 7))
axes_rtlnv = axes_rtlnv.ravel()
coverages_inj = sorted(alkane_inj_com['Target Fractional Surface Coverage'].unique())

for ax, cov, fit_row in zip(axes_rtlnv, coverages_inj, fits_inj_com.itertuples()):
    subset_inj = alkane_inj_com[alkane_inj_com['Target Fractional Surface Coverage'] == cov].dropna(
        subset=['Carbon Number', 'RTlnVg']
    )
    ax.plot(subset_inj['Carbon Number'], subset_inj['RTlnVg'], 'o', markersize=6,
            color='tab:blue', label='Manual (injections)')

    if alkanes is not None:
        subset_sms = alkanes[alkanes['Target Fractional Surface Coverage'] == cov].dropna(
            subset=['Carbon Number', 'RTlnVg']
        )
        ax.plot(subset_sms['Carbon Number'], subset_sms['RTlnVg'], 's', markersize=6,
                color='tab:orange', alpha=0.7, label='SMS (free energy)')

    if len(subset_inj) >= 2 and not np.isnan(fit_row.slope):
        x = subset_inj['Carbon Number'].to_numpy()
        y_fit = fit_row.slope * x + fit_row.intercept
        ax.plot(x, y_fit, '--', linewidth=2, color='tab:blue')
        ax.text(
            0.05, 0.95, f'R² = {fit_row.r_squared:.4f}',
            transform=ax.transAxes, va='top', ha='left',
            fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )

    ax.set_title(f'Coverage = {cov}', fontsize=11)
    ax.set_xlabel('Carbon Number', fontsize=10)
    ax.set_ylabel('RT·ln(V) [J/mol]', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc='lower right')

plt.tight_layout()
plt.savefig(output_dir / 'manual_rtlnv_vs_carbon.png', dpi=150)
print(f"  Saved: {output_dir / 'manual_rtlnv_vs_carbon.png'}")
plt.close(fig_rtlnv)

# 5. Validation against SMS DnG output
print("\n5. VALIDATION VS SMS (DnG COLUMNS)")
print("-" * 60)

if result.dispersive_surface_energy is not None:
    cmp_com = validate_against_sms(gd_com, result.dispersive_surface_energy, retention_type="COM")
    cmp_max = validate_against_sms(gd_max, result.dispersive_surface_energy, retention_type="MAX")

    print("  COM comparison (manual vs SMS DnG & Com):")
    print(cmp_com.to_string(index=False))
    print()

    flagged = cmp_com[cmp_com["exceeds_threshold"]]
    if flagged.empty:
        print("  All coverages within 5% of SMS output.")
    else:
        print(f"  {len(flagged)} coverage(s) exceed 5% deviation:")
        print(flagged[["n/nm", "pct_diff"]].to_string(index=False))

    print("\n  MAX comparison (manual vs SMS DnG & Max):")
    print(cmp_max.to_string(index=False))
else:
    print("  Skipped — dispersive_surface_energy table not present in this CSV.")
    cmp_com = None
    cmp_max = None

# 6. Acid-Base Analysis
print("\n6. ACID-BASE ANALYSIS")
print("-" * 60)

if result.free_energy is not None:
    acid_base = calculate_acid_base_params(result)
    print("Acid-Base Parameters:")
    print(acid_base.to_dataframe().to_string(index=False))
else:
    print("  Skipped — free_energy table required for acid-base analysis.")
    acid_base = None

# 7. Complete Surface Energy Profile
print("\n7. COMPLETE SURFACE ENERGY PROFILE")
print("-" * 60)

if result.free_energy is not None:
    profile = calculate_surface_energy_profile(result)
    print("Surface Energy Profile:")
    print(profile.to_dataframe().to_string(index=False))
    print(f"\n  Exponential Fit Parameters (y = c + a·exp(-b·x)):")
    print(f"  YD:  {profile.fit_params['yd']}")
    print(f"  YAB: {profile.fit_params['yab']}")
    print(f"  YT:  {profile.fit_params['yt']}")
else:
    print("  Skipped — free_energy table required for full surface energy profile.")
    profile = None

# 8. Verification
print("\n8. VERIFICATION")
print("-" * 60)

if profile is not None:
    difference = profile.yt - (profile.yd + profile.yab)
    print(f"YT - (YD + YAB) max difference: {abs(difference).max():.2e}")
    print("Math checks out!" if abs(difference).max() < 1e-10 else "WARNING: YT != YD + YAB")

# 9. Plots
print("\n9. CREATING PLOTS")
print("-" * 60)

# Plot 1: RTlnV vs Carbon Number (SMS pipeline)
if alkanes is not None and fits is not None:
    fig1, axes = plt.subplots(2, 3, figsize=(12, 7))
    axes = axes.ravel()
    coverages = sorted(alkanes['Target Fractional Surface Coverage'].unique())

    for ax, cov, fit_row in zip(axes, coverages, fits.itertuples()):
        subset = alkanes[alkanes['Target Fractional Surface Coverage'] == cov].dropna(
            subset=['Carbon Number', 'RTlnVg']
        )
        ax.plot(subset['Carbon Number'], subset['RTlnVg'], 'o', markersize=6, label='Data')

        if len(subset) >= 2 and not np.isnan(fit_row.slope):
            x = subset['Carbon Number'].to_numpy()
            y_fit = fit_row.slope * x + fit_row.intercept
            ax.plot(x, y_fit, '--', linewidth=2, color='red', label='Linear fit')
            ax.text(
                0.05, 0.95, f'R² = {fit_row.r_squared:.4f}',
                transform=ax.transAxes, va='top', ha='left',
                fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
            )

        ax.set_title(f'Coverage = {cov}', fontsize=11)
        ax.set_xlabel('Carbon Number', fontsize=10)
        ax.set_ylabel('RT·ln(V) [J/mol]', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc='lower right')

    plt.tight_layout()
    plt.savefig(output_dir / 'dorris_gray_sms.png', dpi=150)
    print(f"Saved: {output_dir / 'dorris_gray_sms.png'}")
    plt.close(fig1)

# Plot 2: Manual vs SMS gamma_d comparison
if cmp_com is not None:
    fig2, ax = plt.subplots(figsize=(9, 5))

    ax.plot(cmp_com['n/nm'], cmp_com['gamma_d_manual'], 'o-', label='Manual (COM)', color='tab:blue')
    ax.plot(cmp_com['n/nm'], cmp_com['gamma_d_sms'], 's--', label='SMS DnG & Com', color='tab:orange')

    if cmp_max is not None:
        ax.plot(cmp_max['n/nm'], cmp_max['gamma_d_manual'], '^-', label='Manual (MAX)', color='tab:green')
        ax.plot(cmp_max['n/nm'], cmp_max['gamma_d_sms'], 'D--', label='SMS DnG & Max', color='tab:red')

    ax.set_xlabel('Surface Coverage (n/nm)', fontsize=12)
    ax.set_ylabel('Dispersive Surface Energy (mJ/m²)', fontsize=12)
    ax.set_title('Manual vs SMS Dorris-Gray: γ_s^d vs Coverage', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'manual_vs_sms_gamma_d.png', dpi=150)
    print(f"Saved: {output_dir / 'manual_vs_sms_gamma_d.png'}")
    plt.close(fig2)

# Plot 3: % deviation manual vs SMS
if cmp_com is not None:
    fig3, ax = plt.subplots(figsize=(8, 4))

    ax.bar(
        cmp_com['n/nm'].astype(str),
        cmp_com['pct_diff'],
        label='COM',
        color=['tab:red' if v else 'tab:blue' for v in cmp_com['exceeds_threshold']],
        alpha=0.7,
    )
    ax.axhline(5, color='red', linestyle='--', linewidth=1, label='+5% threshold')
    ax.axhline(-5, color='red', linestyle='--', linewidth=1, label='-5% threshold')
    ax.axhline(0, color='black', linewidth=0.8)

    ax.set_xlabel('Coverage (n/nm)', fontsize=11)
    ax.set_ylabel('% Deviation (manual − SMS) / SMS', fontsize=11)
    ax.set_title('Manual vs SMS Deviation per Coverage Point', fontsize=12)
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'manual_vs_sms_deviation.png', dpi=150)
    print(f"Saved: {output_dir / 'manual_vs_sms_deviation.png'}")
    plt.close(fig3)

# Plot 4: Complete surface energy profile with exponential fits
if profile is not None:
    fig4, ax = plt.subplots(figsize=(10, 6))

    ax.scatter(profile.coverage, profile.yd, label='YD (dispersive)', s=50, color='tab:blue')
    ax.scatter(profile.coverage, profile.yab, label='YAB (acid-base)', s=50, color='tab:orange')
    ax.scatter(profile.coverage, profile.yt, label='YT (total)', s=50, color='tab:green')

    x_fit = np.linspace(0, profile.coverage.max() * 2.0, 200)
    for key, colour in [('yd', 'tab:blue'), ('yab', 'tab:orange'), ('yt', 'tab:green')]:
        p = profile.fit_params[key]
        ax.plot(x_fit, p['c'] + p['a'] * np.exp(-p['b'] * x_fit),
                color=colour, linewidth=2, alpha=0.7)

    ax.set_xlim(0, profile.coverage.max() * 2.0)
    ax.set_xlabel('Surface Coverage (n/nm)', fontsize=12)
    ax.set_ylabel('Surface Energy (mJ/m²)', fontsize=12)
    ax.set_title('Surface Energy Components vs Coverage', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'surface_energy_profile.png', dpi=150)
    print(f"Saved: {output_dir / 'surface_energy_profile.png'}")
    plt.close(fig4)

print("\n" + "=" * 60)
print("ALL TESTS COMPLETE!")
print("=" * 60)
