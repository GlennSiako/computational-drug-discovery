#!/usr/bin/env python3
"""
analyze.py
Pharmacophore analysis using MDAnalysis.

Reads the MD trajectory and detects protein-ligand interactions
at each frame using periodic-boundary-corrected distances.

Interaction types detected:
  Hydrophobic   -- nonpolar C-C contacts within 4.5 A
  Aromatic      -- ring centroid distance within 5.5 A (PHE/TYR/TRP/HIS)
  Hydrogen bond -- N/O donor-acceptor within 3.5 A
  Ionic         -- charged sidechain within 4.0 A of ligand center

All distance calculations use box dimensions so contacts across
periodic boundaries are correctly identified.

Outputs:
  pharmacophore_interactions.csv  -- per-residue frequency table
  pharmacophore_summary.txt       -- plain text summary
"""

import sys
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from collections import Counter
import pandas as pd


print()
print("--- PHARMACOPHORE ANALYSIS ---")
print()

# Load trajectory
if os.path.exists("nvt.gro") and os.path.exists("nvt.xtc"):
    print("  -> Loading MD trajectory (nvt.gro + nvt.xtc)")
    u = mda.Universe("nvt.gro", "nvt.xtc")
    source = "MD trajectory"
elif os.path.exists("em.gro"):
    print("  -> No trajectory found, using minimized structure (em.gro)")
    u = mda.Universe("em.gro")
    source = "Energy-minimized structure"
else:
    print("  No structure found. Run setup.py first.")
    sys.exit(1)

n_frames = len(u.trajectory)
print(f"  Loaded: {len(u.atoms):,} atoms, {n_frames} frames [{source}]")

# Selections
protein = u.select_atoms("protein")

ligand = None
for resname in ["LIG", "BNZ", "BEN", "MOL"]:
    sel = u.select_atoms(f"resname {resname}")
    if len(sel) > 0:
        ligand = sel
        print(f"  Ligand: resname {resname}, {len(sel)} atoms")
        break

if ligand is None:
    ligand = u.select_atoms("not protein and not resname SOL WAT HOH NA CL")
    if len(ligand) == 0:
        print("  Could not find ligand atoms.")
        sys.exit(1)

print(f"  Protein: {len(protein.residues)} residues")


# --- Interaction detectors ---
# All use box=protein.dimensions to apply minimum image convention,
# which correctly handles contacts across periodic boundaries.

def get_hydrophobic(prot, lig, cutoff=4.5):
    """C-C contacts between nonpolar atoms, PBC-corrected."""
    prot_C = prot.select_atoms("name C* and not backbone")
    lig_C  = lig.select_atoms("name C*")
    if len(prot_C) == 0 or len(lig_C) == 0:
        return []
    d = distances.distance_array(
        lig_C.positions, prot_C.positions, box=prot.dimensions
    )
    seen, contacts = set(), []
    for r, c in zip(*np.where(d < cutoff)):
        key = (prot_C[c].resname, prot_C[c].resid)
        if key not in seen:
            seen.add(key)
            contacts.append(f"{prot_C[c].resname}{prot_C[c].resid}")
    return contacts


def get_aromatic(prot, lig, cutoff=5.5):
    """Ring centroid distance, PBC-corrected via minimum image."""
    aromatic_res = {"PHE", "TYR", "TRP", "HIS"}
    lig_C = lig.select_atoms("name C*")
    if len(lig_C) == 0:
        return []
    lig_center = lig_C.positions.mean(axis=0)
    box = prot.dimensions[:3]
    contacts = []
    for res in prot.residues:
        if res.resname in aromatic_res:
            ring = res.atoms.select_atoms("name CG CD* CE* CZ* NE* ND*")
            if len(ring) >= 5:
                diff = ring.positions.mean(axis=0) - lig_center
                diff -= box * np.round(diff / box)
                if np.linalg.norm(diff) < cutoff:
                    contacts.append(f"{res.resname}{res.resid}")
    return contacts


def get_hbonds(prot, lig, cutoff=3.5):
    """N/O donor-acceptor distance, PBC-corrected."""
    lig_acc  = lig.select_atoms("name O* N*")
    prot_don = prot.select_atoms("name N* O*")
    contacts, seen = [], set()
    if len(lig_acc) > 0 and len(prot_don) > 0:
        d = distances.distance_array(
            lig_acc.positions, prot_don.positions, box=prot.dimensions
        )
        for r, c in zip(*np.where(d < cutoff)):
            key = (prot_don[c].resname, prot_don[c].resid)
            if key not in seen:
                seen.add(key)
                contacts.append(f"{prot_don[c].resname}{prot_don[c].resid}(HB)")
    return contacts


def get_ionic(prot, lig, cutoff=4.0):
    """Charged sidechain near ligand center, PBC-corrected."""
    pos_res = {"ARG", "LYS", "HIS"}
    neg_res = {"ASP", "GLU"}
    lig_center = lig.positions.mean(axis=0)
    box = prot.dimensions[:3]
    contacts = []
    for res in prot.residues:
        if res.resname in pos_res | neg_res:
            charged = res.atoms.select_atoms("name NZ NH* NE* OD* OE*")
            if len(charged) > 0:
                diff = charged.positions.mean(axis=0) - lig_center
                diff -= box * np.round(diff / box)
                if np.linalg.norm(diff) < cutoff:
                    sign = "+" if res.resname in pos_res else "-"
                    contacts.append(f"{res.resname}{res.resid}({sign})")
    return contacts


# --- Scan all frames ---
print()
print(f"  -> Scanning {n_frames} frames for interactions...")

all_h, all_ar, all_hb, all_io = [], [], [], []
for ts in u.trajectory:
    all_h.append(get_hydrophobic(protein, ligand))
    all_ar.append(get_aromatic(protein, ligand))
    all_hb.append(get_hbonds(protein, ligand))
    all_io.append(get_ionic(protein, ligand))


# --- Build frequency table ---
def freq_rows(lists, itype):
    counter = Counter(c for frame in lists for c in frame)
    return [
        {
            "residue":          res,
            "interaction_type": itype,
            "frames_detected":  cnt,
            "frequency_%":      round(cnt / n_frames * 100, 1),
        }
        for res, cnt in counter.most_common()
    ]

rows = (
    freq_rows(all_h,  "Hydrophobic") +
    freq_rows(all_ar, "Aromatic (pi-pi)") +
    freq_rows(all_hb, "Hydrogen Bond") +
    freq_rows(all_io, "Ionic")
)

df = pd.DataFrame(rows) if rows else pd.DataFrame(
    columns=["residue", "interaction_type", "frames_detected", "frequency_%"]
)
df = df.sort_values("frequency_%", ascending=False)
df.to_csv("pharmacophore_interactions.csv", index=False)


# --- Print results ---
print()
print("  PHARMACOPHORE FEATURES DETECTED")
print("  " + "-" * 55)
print(f"  {'Residue':<18}  {'Interaction':<18}  Frequency")
print("  " + "-" * 55)

if df.empty:
    print("  No interactions detected.")
else:
    for _, row in df.iterrows():
        res  = str(row["residue"])[:18].ljust(18)
        typ  = str(row["interaction_type"])[:18].ljust(18)
        freq = f"{row['frequency_%']:5.1f}%"
        print(f"  {res}  {typ}  {freq}")

print("  " + "-" * 55)


# --- Write summary file ---
summary = [
    "PHARMACOPHORE SUMMARY",
    "=" * 50,
    f"Source:   {source}",
    f"Frames:   {n_frames}",
    f"Contacts: {len(df)} unique residue contacts",
    "",
]

for itype in ["Hydrophobic", "Aromatic (pi-pi)", "Hydrogen Bond", "Ionic"]:
    sub = df[df["interaction_type"] == itype]
    if not sub.empty:
        summary.append(f"{itype}:")
        for _, row in sub.iterrows():
            summary.append(f"  {row['residue']:<18} {row['frequency_%']:5.1f}%")
        summary.append("")

summary += [
    "=" * 50,
    "DRUG DESIGN INSIGHT:",
    "Benzene binds mainly through HYDROPHOBIC and AROMATIC contacts.",
    "Benzene has no O or N atoms, so hydrogen bond and ionic contacts",
    "are not expected. To design a stronger binder, add groups that",
    "create new interaction types:",
    "  -OH or -NH2   ->  hydrogen bond donor/acceptor",
    "  -COO-         ->  ionic contact with Arg/Lys residues",
    "  Larger ring   ->  more hydrophobic surface area in pocket",
]

with open("pharmacophore_summary.txt", "w") as f:
    f.write("\n".join(summary) + "\n")

print()
print("  Output files written:")
print("    pharmacophore_interactions.csv")
print("    pharmacophore_summary.txt")
print()
print("  Analysis complete.")
print()
