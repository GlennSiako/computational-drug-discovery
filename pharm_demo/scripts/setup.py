#!/usr/bin/env python3
"""
setup.py
Build the complete protein-ligand simulation system.

Steps:
  1. Download 181L.pdb (T4 Lysozyme L99A + benzene)
  2. Generate benzene topology with crystallographic coordinates
  3. Process protein with AMBER99SB-ILDN force field (pdb2gmx)
  4. Patch topol.top with GAFF2 atomtypes and LIG include
  5. Combine protein + ligand coordinates
  6. Create simulation box (dodecahedron)
  7. Solvate with TIP3P water
  8. Add Na+/Cl- ions to neutralize charge
  9. Energy minimization -> em.gro
"""

import subprocess
import sys
import os


def run(cmd, desc, input_text=None):
    print(f"  -> {desc}")
    result = subprocess.run(
        cmd, shell=True, input=input_text,
        text=True, capture_output=True
    )
    if result.returncode != 0:
        print(f"\nFAILED: {desc}")
        print(result.stdout[-3000:])
        print(result.stderr[-3000:])
        sys.exit(1)
    return result.stdout + result.stderr


# Step 1: Download
print()
print("--- STEP 1: Download ---")

if not os.path.exists("181L.pdb"):
    run("curl -fsSL https://files.rcsb.org/download/181L.pdb -o 181L.pdb",
        "Downloading 181L.pdb from RCSB")
    if not os.path.exists("181L.pdb") or os.path.getsize("181L.pdb") < 10000:
        print("  Download failed. Check your internet connection.")
        sys.exit(1)
else:
    print("  -> 181L.pdb already present, skipping download")

with open("181L.pdb") as f:
    lines = f.readlines()

protein_lines = [l for l in lines if l.startswith("ATOM")]
with open("protein_raw.pdb", "w") as f:
    f.writelines(protein_lines)
    f.write("END\n")

bnz_count = sum(1 for l in lines if "BNZ" in l and l.startswith("HETATM"))
print(f"  Protein: {len(protein_lines)} atoms")
print(f"  Benzene in crystal structure: {bnz_count} atoms")


# Step 2: Ligand topology
print()
print("--- STEP 2: Ligand Topology ---")
run("python3 scripts/make_ligand.py", "Generating GAFF2 topology for benzene")


# Step 3: Process protein
print()
print("--- STEP 3: Process Protein ---")
run(
    'echo "6\n1" | gmx pdb2gmx '
    '-f protein_raw.pdb -o protein.gro '
    '-water tip3p -ff amber99sb-ildn -ignh -nobackup 2>&1',
    "pdb2gmx with AMBER99SB-ILDN force field"
)


# Step 4: Patch topology
print()
print("--- STEP 4: Patch Topology ---")

with open("topol.top") as f:
    top = f.read()

gaff2_block = """; GAFF2 atomtypes for benzene (ca = aromatic C, ha = aromatic H)
; This block must appear after forcefield.itp and before any moleculetype.
[ atomtypes ]
; name  at.num   mass      charge  ptype   sigma(nm)      epsilon(kJ/mol)
  ca    6       12.01100  0.000   A       3.39967e-01    3.59824e-01
  ha    1        1.00800  0.000   A       2.59964e-01    6.27600e-02

#include "LIG.itp"

"""

FF_LINE = '#include "amber99sb-ildn.ff/forcefield.itp"'

if "[ atomtypes ]" not in top:
    top = top.replace(
        FF_LINE + "\n\n[ moleculetype ]",
        FF_LINE + "\n\n" + gaff2_block + "[ moleculetype ]"
    )
    top = top.replace(
        "Protein_chain_A     1\n",
        "Protein_chain_A     1\nLIG                 1\n"
    )
    with open("topol.top", "w") as f:
        f.write(top)
    print("  -> GAFF2 atomtypes and LIG.itp inserted into topol.top")
else:
    print("  -> topol.top already patched, skipping")


# Step 5: Combine coordinates
print()
print("--- STEP 5: Combine Coordinates ---")

with open("protein.gro") as f:
    p = f.readlines()
with open("LIG.gro") as f:
    l = f.readlines()

pn = int(p[1].strip())
ln = int(l[1].strip())
combined  = [f"T4 Lysozyme L99A + Benzene\n", f"{pn + ln:5d}\n"]
combined += p[2:-1]
combined += l[2:-1]
combined += [p[-1]]

with open("complex.gro", "w") as f:
    f.writelines(combined)

print(f"  complex.gro: {pn + ln} atoms total")


# Step 6: Simulation box
print()
print("--- STEP 6: Simulation Box ---")
run(
    "gmx editconf -f complex.gro -o box.gro -bt dodecahedron -d 1.0 -nobackup 2>&1",
    "Creating dodecahedron box with 1.0 nm buffer"
)


# Step 7: Solvate
print()
print("--- STEP 7: Solvate ---")
run(
    "gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top -nobackup 2>&1",
    "Adding TIP3P water"
)


# Step 8: Ions
print()
print("--- STEP 8: Add Ions ---")
run(
    "gmx grompp -f config/ions.mdp -c solvated.gro -p topol.top "
    "-o ions.tpr -maxwarn 5 -nobackup 2>&1",
    "grompp for ion placement"
)
run(
    'echo "SOL" | gmx genion -s ions.tpr -o system.gro -p topol.top '
    '-pname NA -nname CL -neutral -nobackup 2>&1',
    "Adding Na+/Cl- ions to neutralize"
)


# Step 9: Energy minimization
print()
print("--- STEP 9: Energy Minimization ---")
run(
    "gmx grompp -f config/em.mdp -c system.gro -p topol.top "
    "-o em.tpr -maxwarn 5 -nobackup 2>&1",
    "grompp for energy minimization"
)
run(
    "gmx mdrun -v -deffnm em -ntmpi 1 -nobackup 2>&1",
    "Running energy minimization"
)

if not os.path.exists("em.gro"):
    print("  em.gro was not created. Check em.log for errors.")
    sys.exit(1)

print(f"  em.gro ready ({os.path.getsize('em.gro'):,} bytes)")
print()
print("  Setup complete. Run: python3 scripts/simulate.py")
print()
