from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


# Resolve paths relative to this file, not the current working directory
HERE = Path(__file__).resolve().parent
CSV_PATH = HERE / "Data" / "Sugar SE.csv"


def _read_rows(path: Path) -> pd.DataFrame:
    # Read raw rows; keep empty strings as empty, don't auto-convert to NaN yet
    return pd.read_csv(path, encoding="cp1252", header=None, dtype=str, keep_default_na=False)


def _is_blank_row(row: pd.Series) -> bool:
    return all((str(v).strip() == "" for v in row.tolist()))


def _maybe_to_numeric(s: pd.Series, *, threshold: float = 0.8) -> pd.Series:
    """
    Convert to numeric only if most non-empty values look numeric.
    Keeps text columns like 'Solvent Name' intact.
    """
    s2 = s.replace("", pd.NA)
    non_empty = s2.dropna()
    if non_empty.empty:
        return s2

    numeric = pd.to_numeric(non_empty, errors="coerce")
    ratio = float(numeric.notna().mean())
    if ratio >= threshold:
        out = pd.to_numeric(s2, errors="coerce")
        return out
    return s2


def extract_section_table(
    rows: pd.DataFrame,
    section_title: str,
    *,
    stop_titles: Optional[List[str]] = None,
    blank_run_to_stop: int = 2,  # require 2 consecutive blank rows to stop
) -> pd.DataFrame:
    stop_titles = stop_titles or []
    col0 = rows[0].astype(str).str.strip()

    # locate section title row
    matches = col0[col0 == section_title]
    if matches.empty:
        raise ValueError(f"Section title not found: {section_title!r}")
    start_idx = int(matches.index[0])

    # header row = next non-blank row after title
    header_idx = None
    for i in range(start_idx + 1, len(rows)):
        if not _is_blank_row(rows.loc[i]):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"No header row found after section {section_title!r}")

    header = [str(x).strip() for x in rows.loc[header_idx].tolist()]
    while header and header[-1] == "":
        header.pop()

    data_start = header_idx + 1

    # find end: stop at next title OR a run of blank rows
    end_idx = len(rows)
    blank_run = 0
    for i in range(data_start, len(rows)):
        if str(rows.loc[i, 0]).strip() in stop_titles:
            end_idx = i
            break

        if _is_blank_row(rows.loc[i]):
            blank_run += 1
            if blank_run >= blank_run_to_stop:
                end_idx = i - blank_run + 1  # end at first blank in the run
                break
        else:
            blank_run = 0

    data = rows.iloc[data_start:end_idx, : len(header)].copy()
    data.columns = header

    # drop fully empty rows/cols
    data = data.replace("", pd.NA)
    data = data.dropna(axis=0, how="all").dropna(axis=1, how="all")

    # numeric conversion
    for c in data.columns:
        data[c] = _maybe_to_numeric(data[c])

    return data


@dataclass(frozen=True)
class IGCResult:
    free_energy: pd.DataFrame
    dispersive_surface_energy: pd.DataFrame
    injection_items: pd.DataFrame
    source_path: Path

    def as_dict(self) -> Dict[str, pd.DataFrame]:
        return {
            "free_energy": self.free_energy,
            "dispersive_surface_energy": self.dispersive_surface_energy,
            "injection_items": self.injection_items,
        }

    def to_csv_dir(self, out_dir: Path) -> None:
        out_dir.mkdir(parents=True, exist_ok=True)
        for name, df in self.as_dict().items():
            df.to_csv(out_dir / f"{name}.csv", index=False)


def parse_sugar_se(path: Path = CSV_PATH) -> IGCResult:
    rows = _read_rows(path)

    free_energy = extract_section_table(
        rows,
        "Free Energy",
        stop_titles=["Dispersive Surface Energy", "Injection Items"],
    )
    dispersive_surface_energy = extract_section_table(
        rows,
        "Dispersive Surface Energy",
        stop_titles=["Injection Items"],
    )
    injection_items = extract_section_table(
        rows,
        "Injection Items",
        stop_titles=[],
    )

    return IGCResult(
        free_energy=free_energy,
        dispersive_surface_energy=dispersive_surface_energy,
        injection_items=injection_items,
        source_path=path,
    )


if __name__ == "__main__":
    result = parse_sugar_se()
    for name, t in result.as_dict().items():
        print(f"{name} {t.shape}")
        print(t.head(3).to_string(index=False))
        print()