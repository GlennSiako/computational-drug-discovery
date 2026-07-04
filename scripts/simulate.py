#!/usr/bin/env python3
"""
simulate.py
Run NVT molecular dynamics simulation at 300 K.

Requires em.gro from setup.py.
Produces nvt.xtc (trajectory) and nvt.gro (final structure).

Simulation parameters:
  Timestep : 0.002 ps
  Steps    : 5000
  Total    : 10 ps
  Temp     : 300 K (V-rescale thermostat)
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
        print(result.stdout[-2000:])
        print(result.stderr[-2000:])
        sys.exit(1)
    return result.stdout + result.stderr


print()
print("--- MD SIMULATION: NVT 300 K, 10 ps ---")
print()

if not os.path.exists("em.gro"):
    print("  em.gro not found. Run setup.py first.")
    sys.exit(1)

run(
    "gmx grompp -f config/nvt.mdp -c em.gro -r em.gro "
    "-p topol.top -o nvt.tpr -maxwarn 5 -nobackup 2>&1",
    "grompp for NVT simulation"
)

print("  -> Running simulation (5000 steps x 0.002 ps = 10 ps)...")
result = subprocess.run(
    "gmx mdrun -v -deffnm nvt -ntmpi 1 -nobackup 2>&1",
    shell=True, text=True, capture_output=True
)

# Print progress lines
for line in (result.stdout + result.stderr).split("\n"):
    if "Step=" in line or "Finished" in line:
        print(f"    {line.strip()}")

if not os.path.exists("nvt.gro"):
    print("  nvt.gro was not created. Check nvt.log for errors.")
    sys.exit(1)

print()
print(f"  nvt.gro  : {os.path.getsize('nvt.gro'):,} bytes")
print(f"  nvt.xtc  : {os.path.getsize('nvt.xtc'):,} bytes")
print()
print("  Simulation complete. Run: python3 scripts/analyze.py")
print()
