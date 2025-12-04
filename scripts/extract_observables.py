#!/usr/bin/env python3
import os
import glob
import numpy as np
import pandas as pd
import mdtraj as md

# --------- CONFIGURATION ---------
TRAJ_NAME = "md_molcompact-nojump-progressive.xtc"
TOP_NAME  = "md_1ns.gro"

# Rg: exclude residues beyond ALA227
# MDTraj resid index is 0-based => ALA227 corresponds to resid <= 226
RG_RESID_CUTOFF = 226


def slugify(s: str) -> str:
    """Make a filesystem-friendly version of a condition name."""
    return (
        s.strip()
         .replace(" ", "_")
         .replace("=", "")
         .replace(",", "")
         .replace(":", "")
    )


def compute_secondary_structure_fractions(traj_protein: md.Trajectory):
    """
    DSSP-based per-frame helix, sheet, and 'other' (coil) fractions.

    Returns:
      helix_frac, sheet_frac, other_frac : 1D arrays of length n_frames
    """
    dssp = md.compute_dssp(traj_protein)  # shape: (n_frames, n_res)
    n_frames, n_res = dssp.shape

    helix_codes = {"H", "G", "I"}
    sheet_codes = {"E", "B"}

    helix_frac = np.zeros(n_frames)
    sheet_frac = np.zeros(n_frames)
    other_frac = np.zeros(n_frames)

    for i in range(n_frames):
        codes = dssp[i]
        h = sum(c in helix_codes for c in codes)
        e = sum(c in sheet_codes for c in codes)
        o = n_res - h - e
        helix_frac[i] = h / n_res
        sheet_frac[i] = e / n_res
        other_frac[i] = o / n_res

    return helix_frac, sheet_frac, other_frac


def process_run(run_dir: str):
    """
    Load trajectory & topology from run_dir and compute per-frame observables.

    Returns:
      pandas.DataFrame with columns:
        time_ns, rg_tail_nm, rmsd_ca_nm, helix_frac, sheet_frac, other_frac
    """
    traj_file = os.path.join(run_dir, TRAJ_NAME)
    top_file  = os.path.join(run_dir, TOP_NAME)

    if not (os.path.isfile(traj_file) and os.path.isfile(top_file)):
        print(f"  [WARN] Missing traj/top in {run_dir}, skipping.")
        return None

    print(f"  Loading {os.path.relpath(run_dir, os.getcwd())} ...")
    traj = md.load(traj_file, top=top_file)

    # Time in ns (mdtraj returns ps)
    time_ns = traj.time / 1000.0

    # Selections
    protein_sel = traj.topology.select("protein")
    ca_sel      = traj.topology.select("protein and name CA")
    rg_sel      = traj.topology.select(f"protein and resid <= {RG_RESID_CUTOFF}")

    traj_protein = traj.atom_slice(protein_sel)
    traj_ca      = traj.atom_slice(ca_sel)
    traj_rg      = traj.atom_slice(rg_sel)

    # Per-frame observables
    rg_tail_nm   = md.compute_rg(traj_rg)             # nm
    rmsd_ca_nm   = md.rmsd(traj_ca, traj_ca[0])       # nm, vs first frame of run
    helix_frac, sheet_frac, other_frac = compute_secondary_structure_fractions(traj_protein)

    df = pd.DataFrame({
        "time_ns":    time_ns,
        "rg_tail_nm": rg_tail_nm,
        "rmsd_ca_nm": rmsd_ca_nm,
        "helix_frac": helix_frac,
        "sheet_frac": sheet_frac,
        "other_frac": other_frac,
    })

    return df


def main():
    # Expect to run from md-protein-stability/data/
    cwd = os.getcwd()
    repo_root = os.path.abspath(os.path.join(cwd, ".."))
    out_root = os.path.join(cwd, "processed")
    os.makedirs(out_root, exist_ok=True)

    # Condition directories: e.g. "300 K 150 mM", "340 K 1000 mM"
    condition_dirs = [d for d in glob.glob("*") if os.path.isdir(d)]
    if not condition_dirs:
        print("[ERROR] No condition directories found in data/.")
        return

    print("Found condition directories:")
    for cond in sorted(condition_dirs):
        print(f"  - {cond}")

    for cond in sorted(condition_dirs):
        cond_path = os.path.join(cwd, cond)
        if not os.path.isdir(cond_path):
            continue

        cond_slug = slugify(cond)
        out_cond_dir = os.path.join(out_root, cond_slug)
        os.makedirs(out_cond_dir, exist_ok=True)

        print(f"\n=== Processing condition: {cond} ===")
        # Run directories: e.g. "r-1", "r-2", "r-3"
        run_dirs = [rd for rd in glob.glob(os.path.join(cond_path, "r-*"))
                    if os.path.isdir(rd)]
        if not run_dirs:
            print(f"  [WARN] No run directories found under {cond}")
            continue

        for run_dir in sorted(run_dirs):
            df = process_run(run_dir)
            if df is None:
                continue

            run_basename = os.path.basename(run_dir)  # e.g. "r-1"
            out_csv = os.path.join(out_cond_dir, f"{run_basename}_observables.csv")

            df.to_csv(out_csv, index=False)
            print(f"  Saved observables to {os.path.relpath(out_csv, repo_root)}")

    print("\nDone. Per-frame observables are in data/processed/<condition_slug>/r-*_observables.csv")


if __name__ == "__main__":
    main()
