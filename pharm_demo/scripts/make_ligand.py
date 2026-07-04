#!/usr/bin/env python3
"""
make_ligand.py
Generate GROMACS topology and coordinates for benzene (GAFF2 force field).

Coordinates are taken directly from the 181L crystal structure so benzene
starts in the correct binding pocket position.

Outputs:
  LIG.gro  -- atom coordinates (nm)
  LIG.itp  -- GAFF2 topology (bonds, angles, dihedrals, charges)

Note: the [ atomtypes ] block for ca/ha is written into topol.top by
setup.py, not here. GROMACS requires atomtypes to appear after the force
field include and before any moleculetype definition.
"""

import math

# Carbon coordinates from 181L.pdb HETATM records (Angstroms)
C_coords = [
    (32.696, -1.329, 18.316),
    (32.130, -0.091, 18.650),
    (30.796,  0.148, 18.241),
    (30.029, -0.840, 17.498),
    (30.595, -2.077, 17.164),
    (31.929, -2.316, 17.573),
]

# Ring center
cx = sum(c[0] for c in C_coords) / 6
cy = sum(c[1] for c in C_coords) / 6
cz = sum(c[2] for c in C_coords) / 6

# Place hydrogens 1.08 A beyond each carbon, radially outward
H_coords = []
for (x, y, z) in C_coords:
    dx, dy, dz = x - cx, y - cy, z - cz
    L = math.sqrt(dx**2 + dy**2 + dz**2)
    s = (L + 1.08) / L
    H_coords.append((cx + dx*s, cy + dy*s, cz + dz*s))

all_coords = C_coords + H_coords
names = [f"C{i+1}" for i in range(6)] + [f"H{i+7}" for i in range(6)]

# --- LIG.gro ---
gro = [
    "Benzene in T4 Lysozyme binding pocket (181L crystal coordinates)",
    "   12",
]
for i, (x, y, z) in enumerate(all_coords):
    # GROMACS uses nm; PDB uses Angstroms
    gro.append(f"    1LIG  {names[i]:>4s}{i+1:5d}{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}")
gro.append("  10.00000  10.00000  10.00000")

with open("LIG.gro", "w") as f:
    f.write("\n".join(gro) + "\n")

# --- LIG.itp ---
CQ, HQ = -0.1306, +0.1306   # AM1-BCC partial charges; sum = 0.0

bonds = [
    (0,1),(1,2),(2,3),(3,4),(4,5),(5,0),    # ring C-C
    (0,6),(1,7),(2,8),(3,9),(4,10),(5,11),  # C-H
]
angles = [
    (5,0,1),(0,1,2),(1,2,3),(2,3,4),(3,4,5),(4,5,0),   # C-C-C
    (6,0,1),(6,0,5),(7,1,0),(7,1,2),(8,2,1),(8,2,3),   # H-C-C
    (9,3,2),(9,3,4),(10,4,3),(10,4,5),(11,5,4),(11,5,0),
]
dihedrals = [
    (0,1,2,3),(1,2,3,4),(2,3,4,5),(3,4,5,0),(4,5,0,1),(5,0,1,2),
    (6,0,1,2),(6,0,1,7),(6,0,5,4),(6,0,5,11),
    (7,1,0,5),(7,1,2,3),(7,1,2,8),
    (8,2,3,4),(8,2,1,0),(8,2,3,9),
    (9,3,4,5),(9,3,2,1),(9,3,4,10),
    (10,4,5,0),(10,4,3,2),(10,4,5,11),
    (11,5,0,1),(11,5,4,3),(11,5,0,6),
]
pairs = [
    (0,2),(0,3),(0,4),(1,3),(1,4),(1,5),(2,4),(2,5),(3,5),
    (0,7),(0,9),(0,10),(1,8),(1,10),(1,11),
    (2,9),(2,11),(2,6),(3,10),(3,6),(3,7),
    (4,11),(4,7),(4,8),(5,6),(5,8),(5,9),
]

itp = [
    "; Benzene (LIG) -- GAFF2 force field",
    "; atomtypes ca and ha are defined in topol.top",
    "",
    "[ moleculetype ]",
    "; name  nrexcl",
    "LIG      3",
    "",
    "[ atoms ]",
    "; nr   type  resnr residue  atom  cgnr     charge     mass",
]
for i in range(6):
    itp.append(f"  {i+1:3d}   ca     1    LIG   {names[i]:4s}    {i+1}   {CQ:.6f}  12.0110")
for i in range(6):
    j = i + 6
    itp.append(f"  {j+1:3d}   ha     1    LIG   {names[j]:4s}    {j+1}   {HQ:.6f}   1.0080")

itp += ["", "[ bonds ]", "; ai  aj  funct  r0(nm)   kb(kJ/mol/nm2)"]
for a, b in bonds:
    params = "0.1398  392459.0" if (a < 6 and b < 6) else "0.1086  307105.0"
    itp.append(f"  {a+1:3d} {b+1:3d}    1   {params}")

itp += ["", "[ angles ]", "; ai  aj  ak  funct  theta0(deg)  ktheta(kJ/mol/rad2)"]
for a, b, c in angles:
    params = "120.000  557.30" if (a < 6 and c < 6) else "120.000  418.40"
    itp.append(f"  {a+1:3d} {b+1:3d} {c+1:3d}    1   {params}")

itp += ["", "[ dihedrals ]", "; ai  aj  ak  al  funct  phase  kd  pn"]
for a, b, c, d in dihedrals:
    itp.append(f"  {a+1:3d} {b+1:3d} {c+1:3d} {d+1:3d}    4   180.000  15.167   2")

itp += ["", "[ pairs ]", "; ai  aj  funct"]
for a, b in pairs:
    itp.append(f"  {a+1:3d} {b+1:3d}    1")

with open("LIG.itp", "w") as f:
    f.write("\n".join(itp) + "\n")

print(f"LIG.gro written  (ligand center: {cx:.1f}, {cy:.1f}, {cz:.1f} A)")
print(f"LIG.itp written  (total charge: {6*CQ + 6*HQ:.6f} e)")
