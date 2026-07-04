#!/bin/bash
set -e

echo ""
echo "Pharmacophore Demo: T4 Lysozyme L99A + Benzene"
echo "================================================"
echo ""

python3 scripts/setup.py    || exit 1
python3 scripts/simulate.py || exit 1
python3 scripts/analyze.py  || exit 1

echo "================================================"
echo "Done."
echo ""
echo "Output files:"
echo "  em.gro                         energy-minimized structure"
echo "  nvt.xtc                        MD trajectory"
echo "  pharmacophore_interactions.csv interaction frequency table"
echo "  pharmacophore_summary.txt      plain text summary"
echo ""
