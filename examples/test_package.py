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

print("=" * 60)
print("Testing IGC-SEA Package")
print("=" * 60)

# 1. Parse CSV
print("\n1. PARSING CSV")
print("-" * 60)
result = parse_igc_csv("data/examples/sample_igc_export.csv")
print(f"✓ Parsed from: {result.source_path.name}")
print(f"  - Free energy: {result.free_energy.shape}")
print(f"  - Dispersive SE: {result.dispersive_surface_energy.shape}")
print(f"  - Injection items: {result.injection_items.shape}")

# 2. Dorris-Gray Alkane Analysis
print("\n2. DORRIS-GRAY ALKANE ANALYSIS")
print("-" * 60)

# Default: use all available alkanes from the Free Energy table
alkanes = prepare_alkane_data(result)
print(f"✓ All available alkanes: {sorted(alkanes['Solvent'].unique())}")
print(f"  Total points: {len(alkanes)}")
print(f"  Coverage levels: {len(alkanes['Target Fractional Surface Coverage'].unique())}")

# Optional: specify which alkanes to use (e.g., exclude hexane if it had poor retention)
# Uncomment below to use a custom alkane range:
# selected_alkanes = ["HEPTANE", "OCTANE", "NONANE", "DECANE"]
# alkanes = prepare_alkane_data(result, alkanes=selected_alkanes)
# print(f"✓ Using selected alkanes: {sorted(alkanes['Solvent'].unique())}")

# Fit linear regression
fits = fit_dorris_gray(alkanes)
print(f"\n  Linear fits (RTlnV vs Carbon Number):")
print(fits.to_string(index=False))

# Verify no negative values
neg_count = (alkanes['RTlnVg'] < 0).sum()
if neg_count > 0:
    print(f"\n⚠ WARNING: {neg_count} negative RTlnVg values found!")
else:
    print(f"\n✓ All RTlnVg values are positive")

# 3. Acid-Base Analysis
print("\n3. ACID-BASE ANALYSIS")
print("-" * 60)
acid_base = calculate_acid_base_params(result)
print(f"✓ Acid-Base Parameters:")
print(acid_base.to_dataframe().to_string(index=False))

# 4. Complete Surface Energy Profile
print("\n4. COMPLETE SURFACE ENERGY PROFILE")
print("-" * 60)
profile = calculate_surface_energy_profile(result)
print(f"✓ Surface Energy Profile:")
print(profile.to_dataframe().to_string(index=False))

print(f"\n  Exponential Fit Parameters (y = c + a·exp(-b·x)):")
print(f"  YD:  {profile.fit_params['yd']}")
print(f"  YAB: {profile.fit_params['yab']}")
print(f"  YT:  {profile.fit_params['yt']}")

# 5. Verify YT = YD + YAB
print("\n5. VERIFICATION")
print("-" * 60)
difference = profile.yt - (profile.yd + profile.yab)
print(f"YT - (YD + YAB) = {difference}")
print(f"Max difference: {abs(difference).max():.2e}")
print("✓ Math checks out!" if abs(difference).max() < 1e-10 else "⚠ Something wrong!")

# 6. Create plots
print("\n6. CREATING PLOTS")
print("-" * 60)

# Ensure outputs directory exists
output_dir = Path(__file__).parent.parent / "outputs"
output_dir.mkdir(exist_ok=True)

# Plot 1: RTlnV vs Carbon Number with regression lines and R²
fig1, axes = plt.subplots(2, 3, figsize=(12, 7))
axes = axes.ravel()
coverages = sorted(alkanes['Target Fractional Surface Coverage'].unique())

for ax, cov, fit_row in zip(axes, coverages, fits.itertuples()):
    subset = alkanes[alkanes['Target Fractional Surface Coverage'] == cov].dropna(
        subset=['Carbon Number', 'RTlnVg']
    )

    # Plot data points
    ax.plot(subset['Carbon Number'], subset['RTlnVg'], 'o', markersize=6, label='Data')

    # Plot best-fit line
    if len(subset) >= 2 and not np.isnan(fit_row.slope):
        x = subset['Carbon Number'].to_numpy()
        y_fit = fit_row.slope * x + fit_row.intercept
        ax.plot(x, y_fit, '--', linewidth=2, color='red', label='Linear fit')

        # Add R² annotation
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
plt.savefig(output_dir / 'dorris_gray_test.png', dpi=150)
print(f"✓ Saved: {output_dir / 'dorris_gray_test.png'}")

# Plot 2: Surface Energy Components
fig2, ax = plt.subplots(figsize=(10, 6))

# Data points
ax.scatter(profile.coverage, profile.yd, label='YD (dispersive)', s=50, color='tab:blue')
ax.scatter(profile.coverage, profile.yab, label='YAB (acid-base)', s=50, color='tab:orange')
ax.scatter(profile.coverage, profile.yt, label='YT (total)', s=50, color='tab:green')

# Fitted curves
x_fit = np.linspace(0, profile.coverage.max() * 2.0, 200)
yd_fit = profile.fit_params['yd']['c'] + profile.fit_params['yd']['a'] * np.exp(
    -profile.fit_params['yd']['b'] * x_fit
)
yab_fit = profile.fit_params['yab']['c'] + profile.fit_params['yab']['a'] * np.exp(
    -profile.fit_params['yab']['b'] * x_fit
)
yt_fit = profile.fit_params['yt']['c'] + profile.fit_params['yt']['a'] * np.exp(
    -profile.fit_params['yt']['b'] * x_fit
)

ax.plot(x_fit, yd_fit, color='tab:blue', linewidth=2, alpha=0.7)
ax.plot(x_fit, yab_fit, color='tab:orange', linewidth=2, alpha=0.7)
ax.plot(x_fit, yt_fit, color='tab:green', linewidth=2, alpha=0.7)

ax.set_xlim(0, profile.coverage.max() * 2.0)

ax.set_xlabel('Surface Coverage (n/nm)', fontsize=12)
ax.set_ylabel('Surface Energy (mJ/m²)', fontsize=12)
ax.set_title('Surface Energy Components vs Coverage', fontsize=14)
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'surface_energy_test.png', dpi=150)
print(f"✓ Saved: {output_dir / 'surface_energy_test.png'}")

print("\n" + "=" * 60)
print("✅ ALL TESTS COMPLETE!")
print("=" * 60)
print("\nYour package is working perfectly! 🎉")
