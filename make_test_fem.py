#!/usr/bin/env python
"""Generate a small but representative GeniE/Sesam FEM test file.

Covers: shells, beams (GBEAMG/GPIPE/GIORH/GBOX/GBARM), truss, solids
(LHEX/TETR/TPRI, one hex deliberately mirrored), spring-to-ground
(GSPR+MGSPRNG), 1-noded mass (GMAS+MGMASS), nodal masses (BNMASS),
fixed + prescribed boundary conditions (BNBCD/BNDISPL), nodal loads
(BNLOAD) and surface loads (BEUSLO on shell and solid face).
"""

def fmt_card(name, fields, cont_lines=()):
    """Format one Sesam FEM card: 8-char name, 16-char fields, 4 per line."""
    out = []
    fields = list(fields)
    first, fields = fields[:4], fields[4:]
    out.append("%-8s" % name + "".join("%16.8E" % f for f in first) + "\n")
    while fields:
        chunk, fields = fields[:4], fields[4:]
        out.append("        " + "".join("%16.8E" % f for f in chunk) + "\n")
    for t in cont_lines:
        out.append("        " + "%-72s" % t + "\n")
    return "".join(out)

cards = []
# --- header ---
cards.append(fmt_card("IDENT", [1, 1, 3, 0]))
# --- nodes: external = internal here ---
coords = {
    # plate corners 1-4 in z=0 plane, beam nodes 5-6 above
    1: (0.0, 0.0, 0.0), 2: (2.0, 0.0, 0.0), 3: (2.0, 1.0, 0.0),
    4: (0.0, 1.0, 0.0), 5: (0.0, 0.0, 3.0), 6: (2.0, 0.0, 3.0),
    # 11-18: LHEX unit cube at x in [10,11] (proper orientation)
    11: (10.0, 0.0, 0.0), 12: (11.0, 0.0, 0.0),
    13: (11.0, 1.0, 0.0), 14: (10.0, 1.0, 0.0),
    15: (10.0, 0.0, 1.0), 16: (11.0, 0.0, 1.0),
    17: (11.0, 1.0, 1.0), 18: (10.0, 1.0, 1.0),
    # 21-28: LHEX unit cube at x in [12,13], mirrored (det J < 0)
    21: (12.0, 0.0, 0.0), 22: (12.0, 1.0, 0.0),
    23: (13.0, 1.0, 0.0), 24: (13.0, 0.0, 0.0),
    25: (12.0, 0.0, 1.0), 26: (12.0, 1.0, 1.0),
    27: (13.0, 1.0, 1.0), 28: (13.0, 0.0, 1.0),
    # 31-34: TETR tetrahedron
    31: (0.0, 10.0, 0.0), 32: (1.0, 10.0, 0.0),
    33: (0.0, 11.0, 0.0), 34: (0.0, 10.0, 1.0),
    # 41-46: TPRI wedge
    41: (0.0, 12.0, 0.0), 42: (1.0, 12.0, 0.0), 43: (0.0, 13.0, 0.0),
    44: (0.0, 12.0, 1.0), 45: (1.0, 12.0, 1.0), 46: (0.0, 13.0, 1.0),
}
for n in sorted(coords):
    cards.append(fmt_card("GNODE", [n, n, 6, 123456]))
for n in sorted(coords):
    x, y, z = coords[n]
    cards.append(fmt_card("GCOORD", [n, x, y, z]))
# --- materials ---
cards.append(fmt_card("MISOSEL", [1, 2.05e11, 0.3, 7850.0, 0, 1.2e-5, 0, 355e6]))
cards.append(fmt_card("MISOSEL", [2, 7.0e10, 0.33, 2700.0, 0, 2.3e-5, 0, 240e6]))
# --- MORSMEL 3: anisotropic membrane material -> MAT2 (Q dropped) ---
# MORSMEL MATNO Q1 Q2 Q3 RHO D11 D21 D22 D31 D32 D33 PS1 PS2 DAMP1 DAMP2 ALPHA1 ALPHA2
cards.append(fmt_card("MORSMEL", [3, 1.0, 0.0, 0.0, 7850.0,
                                  1.62e11, 4.86e10, 1.62e11,
                                  0.0, 0.0, 5.68e10,
                                  0, 0, 0, 0, 1.2e-5, 1.2e-5]))
# --- shell thickness GEONO 100 ---
cards.append(fmt_card("GELTH", [100, 0.012, 0]))
# --- GBEAMG GEONO 200 (I300-ish integral props) ---
cards.append(fmt_card("GBEAMG", [200, 0, 5.35e-3, 8.36e-6, 8.0e-5, 6.2e-6, 0.0,
                                 0, 0, 0, 2.6e-3, 4.3e-3, 0, 0, 0, 0]))
# --- GPIPE GEONO 201 ---
cards.append(fmt_card("GPIPE", [201, 0.46, 0.5, 0.02, 1.0, 1.0, 0, 0]))
# --- GIORH GEONO 202 ---
cards.append(fmt_card("GIORH", [202, 0.3, 0.008, 0.15, 0.012, 0.15, 0.012,
                                 1.0, 1.0, 0, 0, 0]))
# --- GBEAMG GEONO 203 (rod area only) ---
cards.append(fmt_card("GBEAMG", [203, 0, 3.14e-4, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0]))
# --- GBARM GEONO 204 (rectangular bar -> PBEAML BAR) ---
cards.append(fmt_card("GBARM", [204, 0.4, 0.2, 0.2, 1.0, 1.0, 0, 0]))
# --- GBARM GEONO 205 (trapezoidal -> equivalent PBEAM) ---
cards.append(fmt_card("GBARM", [205, 0.4, 0.1, 0.2, 1.0, 1.0, 0, 0]))
# --- GBOX GEONO 206 (equal flanges -> PBEAML BOX) ---
# GBOX GEONO HZ TY TB TT BY SFY SFZ NLOBY NLOBZ
cards.append(fmt_card("GBOX", [206, 0.5, 0.01, 0.02, 0.02, 0.3, 1.0, 1.0, 0, 0]))
# --- GLSEC GEONO 207 (angle -> PBEAML L) ---
# GLSEC GEONO HZ TY BY TZ SFY SFZ K NLOBY NLOBZ
cards.append(fmt_card("GLSEC", [207, 0.25, 0.01, 0.15, 0.012, 1.0, 1.0, 0, 0, 0]))
# --- GCHAN GEONO 208 (channel -> PBEAML CHAN) ---
# GCHAN GEONO HZ TY BY TZ SFY SFZ K NLOBY NLOBZ
cards.append(fmt_card("GCHAN", [208, 0.3, 0.008, 0.12, 0.012, 1.0, 1.0, 0, 0, 0]))
# --- GTONP GEONO 209 (T on plate -> PBEAML I, plate as bottom flange) ---
# GTONP GEONO HZ TY BT TT BP TP SFY SFZ NLOBYT NLOBYB NLOBZ
cards.append(fmt_card("GTONP", [209, 0.35, 0.008, 0.12, 0.012, 0.8, 0.014,
                                1.0, 1.0, 0, 0, 0]))
# --- GUSYI GEONO 210 (unsymmetric I, off-centre web -> PBEAML I + warning) ---
# GUSYI GEONO HZ TY BT B1 TT BB B2 TB SFY SFZ NLOBYT NLOBYB NLOBZ
cards.append(fmt_card("GUSYI", [210, 0.4, 0.01, 0.2, 0.08, 0.015, 0.3, 0.15,
                                0.02, 1.0, 1.0, 0, 0, 0]))
# --- GIORHR GEONO 211 (I with fillets -> PBEAML I, radii dropped) ---
# GIORHR GEONO HZ TY BT TT BB TB SFY SFZ RT RB NLOBYT NLOBYB NLOBZ
cards.append(fmt_card("GIORHR", [211, 0.3, 0.008, 0.15, 0.012, 0.15, 0.012,
                                 1.0, 1.0, 0.02, 0.02, 0, 0, 0]))
# --- GLSECR GEONO 212 (angle with fillet -> PBEAML L, radius dropped) ---
# GLSECR GEONO HZ TY BY TZ SFY SFZ K R NLOBY NLOBZ
cards.append(fmt_card("GLSECR", [212, 0.2, 0.01, 0.12, 0.01,
                                 1.0, 1.0, 0, 0.015, 0, 0]))
# --- GCHANR GEONO 213 (channel with fillet -> PBEAML CHAN, radius dropped) ---
# GCHANR GEONO HZ TY BY TZ SFY SFZ K R NLOBY NLOBZ
cards.append(fmt_card("GCHANR", [213, 0.25, 0.008, 0.1, 0.01,
                                 1.0, 1.0, 0, 0.012, 0, 0]))
# --- GDOBO GEONO 214 (double bottom -> PBEAML I, equal flanges) ---
# GDOBO GEONO HZ TY BY TT TB SFY SFZ NLOBY NLOBZ
cards.append(fmt_card("GDOBO", [214, 1.5, 0.012, 0.9, 0.016, 0.02,
                                1.0, 1.0, 0, 0]))
# --- MGSPRNG 7: diagonal 6x6 + one small off-diagonal K(1,2)=5 ---
# stored row-wise upper triangle: K(1,1..6), K(2,2..6), ..., K(6,6)
cards.append(fmt_card("MGSPRNG", [7, 6,
                                  1.0e6, 5.0, 0, 0, 0, 0,
                                  2.0e6, 0, 0, 0, 0,
                                  3.0e6, 0, 0, 0,
                                  1.0e3, 0, 0,
                                  2.0e3, 0,
                                  3.0e3]))
# --- MGMASS 8: 6x6 diagonal + off-diagonal M(3,1)=7 ---
# stored column-wise lower triangle:
# M(1,1),M(2,1),...,M(6,1), M(2,2),...,M(6,2), ..., M(6,6)
cards.append(fmt_card("MGMASS", [8, 6,
                                 100, 0, 7, 0, 0, 0,
                                 100, 0, 0, 0, 0,
                                 50, 0, 0, 0,
                                 10, 0, 0,
                                 10, 0,
                                 5]))
# --- BNMASS: nodal point masses on node 11 (dofs 1,3,4,6 non-zero) ---
cards.append(fmt_card("BNMASS", [11, 6, 5.0, 0, 5.0, 0.1, 0, 0.1]))
# --- orientation vectors / eccentricity ---
cards.append(fmt_card("GUNIVEC", [300, 0.0, 0.0, 1.0]))
cards.append(fmt_card("GUNIVEC", [301, 0.0, 1.0, 0.0]))
cards.append(fmt_card("GECCEN", [400, 0.1, 0.0, 0.05]))
# --- GECC 401: local-frame eccentricity (IOPT=2: EX EY EZ) ---
cards.append(fmt_card("GECC", [401, 2, 0.0, 0.05, 0.1]))
# --- BELFIX 500: axial release both ends (A1=0, rest fixed) ---
cards.append(fmt_card("BELFIX", [500, 1, 0, 0, 0.0, 1.0, 1.0, 1.0,
                                 1.0, 1.0, 0.0]))
# --- BELFIX 501: fractional fixation on dof 6 -> warning path ---
cards.append(fmt_card("BELFIX", [501, 1, 0, 0, 1.0, 1.0, 1.0, 1.0,
                                 1.0, 0.5, 0.0]))
# --- user text -> BDF comment ---
cards.append(fmt_card("TEXT", [1, 0, 1, 28], ["fem2bdf regression model"]))
# --- elements ---
def elem(elno, eltyp, nodes, matno, geono=0, eccno=0, transno=0, fixno=0):
    cards.append(fmt_card("GELMNT1", [elno, elno, eltyp, 0] + list(nodes)))
    cards.append(fmt_card("GELREF1", [elno, matno, 0, 0, 0, 0, 0, 0,
                                      geono, fixno, eccno, transno]))

# E1 quad (FQUS=24) on thickness 100 / mat 1
elem(1, 24, [1, 2, 3, 4], 1, geono=100)
# E2 tria (FTRS=25) on thickness 100 / mat 2
elem(2, 25, [1, 4, 5], 2, geono=100)
# E3 beam (BEAS=15) nodes 5-6, GBEAMG 200 / mat 1 / trans 300
#   / BELFIX 500 (axial release both ends -> PA=1 PB=1)
elem(3, 15, [5, 6], 1, geono=200, transno=300, fixno=500)
# E4 beam nodes 6-1, GPIPE 201 / mat 2 / trans 301 / ecc 400
elem(4, 15, [6, 1], 2, geono=201, eccno=400, transno=301)
# E5 beam nodes 5-1, GIORH 202 / mat 1 / trans 300
#   / BELFIX 501 (fractional fixation -> warning, no pin flags)
elem(5, 15, [5, 1], 1, geono=202, transno=300, fixno=501)
# E6 truss (TESS=10) nodes 1-6, GBEAMG 203 / mat 1
elem(6, 10, [1, 6], 1, geono=203)
# E7 LHEX (21) proper cube, mat 1
elem(7, 21, [11, 12, 13, 14, 15, 16, 17, 18], 1)
# E8 LHEX (21) mirrored cube, mat 2 (converter must flip node order)
elem(8, 21, [21, 22, 23, 24, 25, 26, 27, 28], 2)
# E9 TETR (33), mat 1
elem(9, 33, [31, 32, 33, 34], 1)
# E10 TPRI (32), mat 2
elem(10, 32, [41, 42, 43, 44, 45, 46], 2)
# E11 GSPR (18) spring to ground at node 5, MGSPRNG 7
elem(11, 18, [5], 7)
# E12 GMAS (11) 1-noded mass at node 6, MGMASS 8
elem(12, 11, [6], 8)
# E13 beam nodes 6-11, GBARM 204 (rectangular -> PBEAML BAR)
#   / GECC 401 (local eccentricity, rotated to global)
elem(13, 15, [6, 11], 1, geono=204, eccno=401, transno=300)
# E14 beam nodes 6-12, GBARM 205 (trapezoidal -> equiv PBEAM)
elem(14, 15, [6, 12], 1, geono=205, transno=300)
# E15 beam nodes 12-16, GBOX 206 (-> PBEAML BOX)
elem(15, 15, [12, 16], 1, geono=206, transno=300)
# E16 beam nodes 13-14, GLSEC 207 (angle -> PBEAML L)
elem(16, 15, [13, 14], 1, geono=207, transno=300)
# E17 beam nodes 14-18, GCHAN 208 (channel -> PBEAML CHAN)
elem(17, 15, [14, 18], 1, geono=208, transno=300)
# E18 beam nodes 16-17, GTONP 209 (T on plate -> PBEAML I)
elem(18, 15, [16, 17], 1, geono=209, transno=300)
# E19 beam nodes 17-18, GUSYI 210 (unsymmetric I, off-centre web)
elem(19, 15, [17, 18], 1, geono=210, transno=300)
# E20 beam nodes 15-18, GIORHR 211 (I with fillets)
elem(20, 15, [15, 18], 1, geono=211, transno=300)
# E21 beam nodes 15-16, GLSECR 212 (angle with fillet)
elem(21, 15, [15, 16], 1, geono=212, transno=300)
# E22 beam nodes 11-15, GCHANR 213 (channel with fillet)
elem(22, 15, [11, 15], 1, geono=213, transno=300)
# E23 beam nodes 22-23, GDOBO 214 (double bottom -> PBEAML I)
elem(23, 15, [22, 23], 1, geono=214, transno=300)
# E24 tria (FTRS=25) nodes 1-2-5, thickness 100 / MORSMEL mat 3
elem(24, 25, [1, 2, 5], 3, geono=100)
# E25 quad nodes 1-2-6-5, GEONO 101 has NO GELTH card
#   -> default PSHELL t=1.0 emitted with warning
elem(25, 24, [1, 2, 6, 5], 1, geono=101)
# E26 beam nodes 2-3, GEONO 215 has NO section card at all
#   -> property skipped, element skipped (no dangling PID)
elem(26, 15, [2, 3], 1, geono=215, transno=300)
# --- set definitions: TDSETNAM + GSETMEMB (informational comments) ---
cards.append(fmt_card("TDSETNAM", [4, 5, 104, 0], ["SetA"]))
cards.append(fmt_card("GSETMEMB", [7, 5, 1, 1, 0, 1, 2]))
cards.append(fmt_card("GSETMEMB", [7, 6, 1, 2, 0, 7, 8]))
# --- linear dependence: node 4 dof 3 slaved to node 3 dof 3 (b=1) ---
# BLDEP NODENO CNOD NDDOF NDEP / DEPDOF INDEPDOF b (padded to 9 fields)
cards.append(fmt_card("BLDEP", [4, 3, 6, 1, 3, 3, 1.0, 0, 0]))
cards.append(fmt_card("BNBCD", [4, 6, 0, 0, 3, 0, 0, 0]))  # node4 dof3 dependent
# --- boundary conditions ---
cards.append(fmt_card("BNBCD", [1, 6, 1, 1, 1, 1, 1, 1]))  # node1: 123456 fixed
cards.append(fmt_card("BNBCD", [6, 6, 1, 1, 1, 0, 0, 0]))  # node6: 123 fixed
# node12: dofs 1,2 prescribed (values on BNDISPL)
cards.append(fmt_card("BNBCD", [12, 6, 2, 2, 0, 0, 0, 0]))
# node13: dof 3 prescribed but no BNDISPL (exercises the warning path)
cards.append(fmt_card("BNBCD", [13, 6, 0, 0, 2, 0, 0, 0]))
# prescribed values for load case 1, node 12: d1=0.01, d2=-0.02
# layout: LLC DTYPE COMPLX <pad> NODENO NDOF RDISP1..6
cards.append(fmt_card("BNDISPL", [1, 1, 0, 0, 12, 6,
                                  0.01, -0.02, 0, 0, 0, 0]))
# --- load case names (NFIELD=4, CODNAM=100+namelen) ---
cards.append(fmt_card("TDLOAD", [4, 1, 103, 0], ["LC1"]))
cards.append(fmt_card("TDLOAD", [4, 2, 103, 0], ["LC2"]))
# --- nodal loads ---
# BNLOAD layout as written by GeniE: first line padded with a spare
# zero field after COMPLX, NODENO starts the second line.
# LLC1: node5 force (1000, 0, -500)
cards.append(fmt_card("BNLOAD", [1, 1, 0, 0, 5, 6, 1000.0, 0.0, -500.0, 0.0, 0.0, 0.0]))
# LLC2: node6 moment (0, 250, 0) + force z
cards.append(fmt_card("BNLOAD", [2, 1, 0, 0, 6, 6, 0.0, 0.0, 100.0, 0.0, 250.0, 0.0]))
# --- surface loads ---
# BEUSLO: LLC LOTYP COMPLX LAYER ELNO NDOF INTNO SIDE RLOAD1..
# LC1: uniform 1000 on top face (SIDE=3 -> -p) of quad E1
cards.append(fmt_card("BEUSLO", [1, 1, 0, 0, 1, 4, 0, 3,
                                 1000.0, 1000.0, 1000.0, 1000.0]))
# LC2: varying pressure on bottom face (SIDE=1 -> +p) of quad E1
cards.append(fmt_card("BEUSLO", [2, 1, 0, 0, 1, 4, 0, 1,
                                 500.0, 600.0, 700.0, 800.0]))
# LC1: uniform 2000 on side 4 (SIDE=000100) of hex E7
cards.append(fmt_card("BEUSLO", [1, 1, 0, 0, 7, 4, 0, 100,
                                 2000.0, 2000.0, 2000.0, 2000.0]))
# LC2: follower pressure (LOTYP=-1) on top face of tria E2 (warning)
cards.append(fmt_card("BEUSLO", [2, -1, 0, 0, 2, 3, 0, 3,
                                 900.0, 900.0, 900.0]))
# LC2: component-form load (LOTYP=2) on E1 -> skipped with warning
cards.append(fmt_card("BEUSLO", [2, 2, 0, 0, 1, 3, 0, 0,
                                 1.0, 2.0, 3.0]))
# LC2: hydro-pressure marker (LOTYP=3) on E7 -> skipped with warning
cards.append(fmt_card("BEUSLO", [2, 3, 0, 0, 7, 1, 0, 0, 1.0]))
# --- end of superelement ---
cards.append(fmt_card("IEND", [0, 0, 0, 0]))

with open("test_model.fem", "w") as f:
    f.write("".join(cards))
print("test_model.fem written,", len(cards), "cards")
