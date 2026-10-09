
#!/usr/bin/env python3
"""
Detect mixed peaks in .ab1 Sanger chromatograms using PLOC-only detection.

Features:
  - Uses ABI-provided peak index (PLOC2) for each base.
  - Primary and secondary heights taken at the exact peak index.
  - Baseline correction: subtract local median background per channel.
  - Quality filter: skip low-quality bases (Phred < min-quality).
  - Primary strength guard: require sufficiently strong primary peak.
  - Per-file CSV + single summary CSV ("none detected" if no mixtures).
  - Optional IUPAC FASTA and chromatogram snapshots.
  - Debug mode for inspecting corrected heights/ratios.
  - Optional override for channel base order.

Dependencies: biopython, numpy, pandas, matplotlib
Install:     python -m pip install biopython numpy pandas matplotlib
"""

import argparse
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

# IUPAC codes for degenerate bases
IUPAC_MAP = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
    frozenset(["A", "C", "G"]): "V",
    frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "C", "G", "T"]): "N",
}


def _to_array(v) -> Optional[np.ndarray]:
    if v is None:
        return None
    try:
        return np.array(v, dtype=np.int32)
    except Exception:
        if isinstance(v, (bytes, bytearray)):
            try:
                return np.frombuffer(v, dtype=">u2")
            except Exception:
                return None


def _get_base_order(abif_raw: Dict, override: Optional[str]) -> str:
    """Return base order (e.g., 'GATC'), using override if provided."""
    if override:
        override = override.strip().upper()
        if len(override) == 4 and set(override) == set("ACGT"):
            return override
        else:
            raise ValueError(f"Invalid --base-order '{override}'. Must be a permutation of ACGT.")

    order = abif_raw.get("FWO_")
    if isinstance(order, (bytes, bytearray)):
        try:
            s = order.decode("ascii", errors="ignore").strip().upper()
            if len(s) == 4 and set(s) == set("ACGT"):
                return s
        except Exception:
            pass
    if isinstance(order, str):
        s = order.strip().upper()
        if len(s) == 4 and set(s) == set("ACGT"):
            return s

    return "GATC"  # fallback


def _find_trace_channels(abif_raw: Dict, base_order_override: Optional[str]) -> Dict[str, np.ndarray]:
    candidates = []
    for k, v in abif_raw.items():
        if isinstance(k, str) and k.startswith("DATA"):
            arr = _to_array(v)
            if arr is not None and arr.ndim == 1 and arr.size > 100:
                candidates.append((k, arr))
    if not candidates:
        raise RuntimeError("No DATA* trace arrays found in ABIF.")
    max_len = max(arr.size for _, arr in candidates)
    channels = [(k, arr) for k, arr in candidates if arr.size == max_len]

    def suffix_num(tag: str) -> int:
        try:
            return int(tag[4:])
        except Exception:
            return -1

    if len(channels) >= 4:
        channels = sorted(channels, key=lambda kv: suffix_num(kv[0]))[-4:]
    if len(channels) != 4:
        raise RuntimeError(f"Found {len(channels)} channels; expected 4.")
    channels_sorted = sorted(channels, key=lambda kv: suffix_num(kv[0]))
    base_order = _get_base_order(abif_raw, base_order_override)
    return {base: arr for base, (_, arr) in zip(base_order, channels_sorted)}


def _get_sample_name(abif_raw: Dict, path: Path) -> str:
    for key in ("SMPL1", "TUBE1", "RunName", "RUND1"):
        val = abif_raw.get(key)
        if isinstance(val, (bytes, bytearray)):
            try:
                s = val.decode("ascii", errors="ignore").strip()
                if s:
                    return s
            except Exception:
                pass
        elif isinstance(val, str) and val.strip():
            return val.strip()
    return path.stem


def local_background(arr: np.ndarray, idx: int, bg_window: int) -> float:
    left = max(0, idx - bg_window)
    right = min(arr.size, idx + bg_window + 1)
    seg = arr[left:right]
    return float(np.median(seg)) if seg.size else 0.0


def corrected_height(arr: np.ndarray, idx: int, bg_window: int) -> float:
    h = arr[idx]
    bg = local_background(arr, idx, bg_window)
    return float(max(0.0, h - bg))


def plot_snapshot(trace_map: Dict[str, np.ndarray], idx: int, window: int, out_path: Path, title: str, called: str, secondary: str):
    left = max(0, idx - window)
    right = min(idx + window + 1, min(len(trace_map["A"]), len(trace_map["C"]), len(trace_map["G"]), len(trace_map["T"])))
    plt.figure(figsize=(9, 4.8))
    colors = {"A": "green", "C": "blue", "G": "black", "T": "red"}
    for base in "ACGT":
        y = trace_map[base][left:right]
        x = np.arange(left, right)
        plt.plot(x, y, label=base, color=colors[base], linewidth=1.2)
    plt.axvline(idx, color="gray", linestyle="--", linewidth=1)
    plt.title(f"{title} | called={called}, secondary={secondary}")
    plt.xlabel("Trace index")
    plt.ylabel("Signal (RFU)")
    plt.legend(loc="upper right", ncol=4, frameon=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close()


def analyze_ab1_file(path: Path, ratio: float, min_secondary_rfu: int, min_primary_rfu: int, bg_window: int, ignore_ends: int, min_quality: int, write_iupac: bool, snapshots: bool, snapshot_window: int, base_order_override: Optional[str], out_dir: Path) -> Dict:
    record = SeqIO.read(str(path), "abi")
    abif = record.annotations.get("abif_raw", {})
    seq_str = str(record.seq)
    ppositions = abif.get("PLOC2") or abif.get("PLOC1")
    if ppositions is None:
        raise RuntimeError("Missing PLOC array.")
    ppositions = list(ppositions)
    qualities = record.letter_annotations.get("phred_quality", [None] * len(seq_str))
    trace_map = _find_trace_channels(abif, base_order_override)
    sample_name = _get_sample_name(abif, path)

    rows = []
    iupac_seq = list(seq_str)
    total_bases = len(seq_str)
    start_idx = ignore_ends
    end_idx = total_bases - ignore_ends
    snapshot_dir = out_dir / f"{path.stem}_snapshots"
    if snapshots:
        snapshot_dir.mkdir(parents=True, exist_ok=True)
    mixed_positions = []

    for i in range(total_bases):
        called = seq_str[i].upper()
        if called not in "ACGT" or i < start_idx or i >= end_idx:
            continue
        q = qualities[i] if qualities and i < len(qualities) else None
        if q is not None and q < min_quality:
            continue
        idx = int(ppositions[i])
        heights_corr = {b: corrected_height(trace_map[b], idx, bg_window) for b in "ACGT"}
        primary_h = heights_corr.get(called, 0.0)
        if primary_h < min_primary_rfu:
            continue
        passing_secondaries = {b: h for b, h in heights_corr.items() if b != called and h >= min_secondary_rfu and (h / primary_h) >= ratio}
        if passing_secondaries:
            sec_base, sec_h = max(passing_secondaries.items(), key=lambda kv: kv[1])
            sec_ratio = sec_h / primary_h
            rows.append({"file": path.name, "sample_name": sample_name, "base_index_1based": i + 1, "called_base": called, "secondary_base": sec_base, "primary_height_corrected": round(primary_h, 2), "secondary_height_corrected": round(sec_h, 2), "secondary_primary_ratio": round(sec_ratio, 3), "peak_location_index": idx, "phred_quality": q})
            mixed_positions.append(str(i + 1))
            degenerate_set = frozenset([called] + list(passing_secondaries.keys()))
            iupac_char = IUPAC_MAP.get(degenerate_set, called)
            if write_iupac:
                iupac_seq[i] = iupac_char
            if snapshots:
                out_path = snapshot_dir / f"pos_{i+1}_{called}_{sec_base}.png"
                plot_snapshot(trace_map, idx, snapshot_window, out_path, f"{path.stem} pos {i+1}", called, sec_base)

    out_dir.mkdir(parents=True, exist_ok=True)
    per_file_csv = out_dir / f"{path.stem}_mixed_positions.csv"
    pd.DataFrame(rows).to_csv(per_file_csv, index=False)
    fasta_path = None
    if write_iupac:
        fasta_path = out_dir / f"{path.stem}_iupac.fasta"
        with open(fasta_path, "w") as fh:
            fh.write(f">{sample_name}_iupac\n")
            fh.write("".join(iupac_seq) + "\n")
    summary = {"file": path.name, "sample_name": sample_name, "total_bases": total_bases, "mixed_count": len(rows), "mixed_pct": round(100.0 * len(rows) / max(total_bases, 1), 2), "mixed_positions": ", ".join(mixed_positions) if mixed_positions else "none detected", "per_file_csv": str(per_file_csv), "iupac_fasta": str(fasta_path) if fasta_path else None}
    return summary


def main():
    parser = argparse.ArgumentParser(description="Detect mixed peaks using PLOC-only detection.")
    parser.add_argument("input", type=str, help="Input folder containing .ab1 files.")
    parser.add_argument("-o", "--out", type=str, default="ab1_mixed_out", help="Output folder.")
    parser.add_argument("-r", "--ratio", type=float, default=0.20, help="Secondary/primary ratio threshold.")
    parser.add_argument("--strict", action="store_true", help="Set ratio=0.50.")
    parser.add_argument("--lenient", action="store_true", help="Set ratio=0.25.")
    parser.add_argument("--min-secondary-rfu", type=int, default=100, help="Minimum secondary RFU.")
    parser.add_argument("--min-primary-rfu", type=int, default=200, help="Minimum primary RFU.")
    parser.add_argument("--bg-window", type=int, default=25, help="Background estimation window.")
    parser.add_argument("--ignore-ends", type=int, default=20, help="Ignore ends.")
    parser.add_argument("--min-quality", type=int, default=20, help="Minimum Phred quality.")
    parser.add_argument("--recursive", action="store_true", help="Search recursively.")
    parser.add_argument("--base-order", type=str, default=None, help="Override base order (e.g., GATC).")
    parser.add_argument("--iupac", action="store_true", help="Write IUPAC FASTA.")
    parser.add_argument("--snapshots", action="store_true", help="Generate chromatogram snapshots.")
    parser.add_argument("--snapshot-window", type=int, default=50, help="Snapshot window.")
    args = parser.parse_args()
    if args.strict:
        args.ratio = 0.50
    if args.lenient:
        args.ratio = 0.25
    in_dir = Path(args.input)
    out_dir = Path(args.out)
    if not in_dir.exists():
        raise SystemExit(f"Input folder not found: {in_dir}")
    pattern_func = Path.rglob if args.recursive else Path.glob
    files = list(pattern_func(in_dir, "*.ab1"))
    if not files:
        raise SystemExit("No .ab1 files found.")
    summaries = []
    for f in files:
        try:
            summary = analyze_ab1_file(f, args.ratio, args.min_secondary_rfu, args.min_primary_rfu, args.bg_window, args.ignore_ends, args.min_quality, args.iupac, args.snapshots, args.snapshot_window, args.base_order, out_dir)
            summaries.append(summary)
            print(f"[OK] {summary['sample_name']} ({f.name}): {summary['mixed_count']} mixed -> {summary['mixed_positions']}")
        except Exception as e:
            print(f"[ERROR] {f.name}: {e}")
    summary_df = pd.DataFrame(summaries)
    summary_csv = out_dir / "summary_mixed_peaks.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"\nSummary written to: {summary_csv}")


if __name__ == "__main__":
    main()
