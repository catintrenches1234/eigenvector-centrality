"""
Molecular Graph Definitions
===========================
Carbon skeleton graphs for alkane isomers (C4-C8).
"""

MOLECULES = {
    # ── Butane isomers (C4H10) ──
    "n-butane": {
        "n": 4,
        "edges": [(0,1),(1,2),(2,3)],
        "formula": "C4H10",
        "note": "Linear chain: C1-C2-C3-C4"
    },
    "isobutane (2-methylpropane)": {
        "n": 4,
        "edges": [(0,1),(1,2),(1,3)],
        "formula": "C4H10",
        "note": "Star shape: C2 is the central carbon (degree 3)"
    },

    # ── Pentane isomers (C5H12) ──
    "n-pentane": {
        "n": 5,
        "edges": [(0,1),(1,2),(2,3),(3,4)],
        "formula": "C5H12",
        "note": "Linear chain"
    },
    "isopentane (2-methylbutane)": {
        "n": 5,
        "edges": [(0,1),(1,2),(2,3),(2,4)],
        "formula": "C5H12",
        "note": "Branch at C3"
    },
    "neopentane (2,2-dimethylpropane)": {
        "n": 5,
        "edges": [(0,1),(1,2),(1,3),(1,4)],
        "formula": "C5H12",
        "note": "Star: C2 connected to 4 carbons"
    },

    # ── Hexane isomers (C6H14) ──
    "n-hexane": {
        "n": 6,
        "edges": [(0,1),(1,2),(2,3),(3,4),(4,5)],
        "formula": "C6H14",
        "note": "Linear chain"
    },
    "2-methylpentane": {
        "n": 6,
        "edges": [(0,1),(1,2),(2,3),(3,4),(1,5)],
        "formula": "C6H14",
        "note": "Methyl branch at C2"
    },
    "3-methylpentane": {
        "n": 6,
        "edges": [(0,1),(1,2),(2,3),(3,4),(2,5)],
        "formula": "C6H14",
        "note": "Methyl branch at C3 (middle)"
    },
    "2,2-dimethylbutane": {
        "n": 6,
        "edges": [(0,1),(1,2),(2,3),(1,4),(1,5)],
        "formula": "C6H14",
        "note": "Two methyl branches at C2 (degree 4 node)"
    },
    "2,3-dimethylbutane": {
        "n": 6,
        "edges": [(0,1),(1,2),(2,3),(1,4),(2,5)],
        "formula": "C6H14",
        "note": "One methyl branch each at C2 and C3"
    },

    # ── n-Octane (C8H18) for reference ──
    "n-octane": {
        "n": 8,
        "edges": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7)],
        "formula": "C8H18",
        "note": "8-carbon linear chain (Rani & Balamurugan reference)"
    },
}

# Known boiling points (°C) for comparison
BOILING_POINTS = {
    "n-butane": -0.5,
    "isobutane (2-methylpropane)": -11.7,
    "n-pentane": 36.1,
    "isopentane (2-methylbutane)": 27.7,
    "neopentane (2,2-dimethylpropane)": 9.5,
    "n-hexane": 68.7,
    "2-methylpentane": 60.3,
    "3-methylpentane": 63.3,
    "2,2-dimethylbutane": 49.7,
    "2,3-dimethylbutane": 58.0,
    "n-octane": 125.7,
}
