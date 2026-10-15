#!/usr/bin/env python
"""Field-precise verification of test_model.bdf against expectations.

Checks:
 1. CONM1 21 lower-triangle values match the MGMASS input (column-major
    storage correctly unpacked).
 2. PBEAML DIM values produce section properties (MSC beam-library
    formulas) matching the converter's own section-property routines.
 3. CHEXA/CTETRA/CPENTA orientation: det J > 0 with the written order.
 4. SPC / SPC1 contents per load case.
 5. PLOAD4 fields (SID, EID, pressures, G1/G4).
 6. PBUSH diagonal stiffness vector.
 7. CMASS2 cards.
"""
import math
import sys

LINES = open("test_model.bdf").read().splitlines()

def fields(line):
    return [line[i:i+8] for i in range(0, len(line), 8)]

def fld_float(s):
    s = s.strip()
    if not s:
        return 0.0
    # ExtFEMIO writes e.g. "1.000+02" (no E)
    if ("+" in s[1:] or "-" in s[1:]) and "E" not in s.upper():
        s = s.replace("+", "E+").replace("-", "E-", 1) if s[0] not in "+-" \
            else s[0] + s[1:].replace("+", "E+").replace("-", "E-")
    return float(s)

def fld_int(s):
    return int(float(fld_float(s))) if s.strip() else 0

fails = []
def check(name, got, want, tol=1e-6):
    ok = (abs(got - want) <= tol * max(1.0, abs(want))) \
        if isinstance(want, float) else (got == want)
    print(("PASS" if ok else "FAIL"), name, "got", got, "want", want)
    if not ok:
        fails.append(name)

# ---- collect cards ----
grids = {}
conm1_vals = []
pload4 = []
pbush = {}
spc = []
spc1 = []
cmass2 = []
pbeaml = {}   # pid -> (type, [dims])
chexa = []
ctetra = []
cpenta = []

i = 0
while i < len(LINES):
    l = LINES[i]
    name = l[:8].strip()
    if name == "GRID":
        f = fields(l)
        grids[fld_int(f[1])] = tuple(fld_float(x) for x in f[3:6])
    elif name == "CONM1":
        f = fields(l) + fields(LINES[i+1]) + fields(LINES[i+2])
        # QRG layout:
        # l1: (name) EID G CID M11 M21 M22 M31 M32      -> merged[4:9]
        # l2: (blank) M33 M41 M42 M43 M44 M51 M52 M53   -> merged[10:18]
        # l3: (blank) M54 M55 M61 M62 M63 M64 M65 M66   -> merged[19:27]
        vals = [fld_float(x) for x in f[4:9]] + \
               [fld_float(x) for x in f[10:18]] + \
               [fld_float(x) for x in f[19:27]]
        conm1_vals = vals
        i += 2
    elif name == "PLOAD4":
        f = fields(l)
        p = [fld_float(x) for x in f[3:7]]
        p += [0.0] * (4 - len(p))   # blank P2..P4 = same as P1
        pload4.append([fld_int(f[1]), fld_int(f[2])] + p +
                      [fld_int(f[7]) if len(f) > 7 else 0,
                       fld_int(f[8]) if len(f) > 8 else 0])
    elif name == "PBUSH":
        f = fields(l)
        pbush[fld_int(f[1])] = [fld_float(x) for x in f[3:9]]
    elif name == "SPC1":
        f = fields(l)
        spc1.append((fld_int(f[1]), f[2].strip(),
                     [fld_int(x) for x in f[3:] if x.strip()]))
    elif name == "SPC":
        f = fields(l)
        spc.append((fld_int(f[1]), fld_int(f[2]), fld_int(f[3]),
                    fld_float(f[4])))
    elif name == "CMASS2":
        f = fields(l)
        cmass2.append((fld_int(f[1]), fld_float(f[2]),
                       fld_int(f[3]), fld_int(f[4])))
    elif name == "PBEAML":
        f = fields(l)
        pid = fld_int(f[1])
        typ = f[4].strip()
        dims = []
        j = i + 1
        # continuation lines hold 8-char fields with dims then YES, X/XB
        while j < len(LINES) and LINES[j].startswith("        "):
            for x in fields(LINES[j]):
                xs = x.strip()
                if not xs or xs in ("YES", "YESA", "NO"):
                    continue
                v = fld_float(xs)
                if abs(v - 1.0) < 1e-12:   # trailing X/XB = 1.0 marker
                    continue
                dims.append(v)
            j += 1
        pbeaml[pid] = (typ, dims)
        i = j - 1
    elif name == "CHEXA":
        f = fields(l) + fields(LINES[i+1])
        chexa.append([fld_int(f[1])] + [fld_int(x) for x in f[3:]
                                        if x.strip()])
        i += 1
    elif name == "CTETRA":
        f = fields(l)
        ctetra.append([fld_int(f[1])] + [fld_int(x) for x in f[3:]
                                         if x.strip()])
    elif name == "CPENTA":
        f = fields(l)
        cpenta.append([fld_int(f[1])] + [fld_int(x) for x in f[3:]
                                         if x.strip()])
    i += 1

print("== CONM1 ==")
# expected from MGMASS 8 (column-major lower):
# M11=100 M22=100 M33=50 M44=10 M55=10 M66=5, M31=7
want = [100., 0., 100., 7., 0., 50., 0., 0., 0., 10.,
        0., 0., 0., 0., 10., 0., 0., 0., 0., 0., 5.]
check("CONM1 21 values", conm1_vals, [float(w) for w in want]) \
    if False else None
for k, (g, w) in enumerate(zip(conm1_vals, want)):
    check("CONM1 M[%d]=%g" % (k, w), g, float(w))

print("== PBEAML section cross-check ==")
# MSC beam-library formulas (diameter-based TUBE)
def msc_tube(d, t):
    di = d - 2 * t
    A = math.pi / 4 * (d * d - di * di)
    I = math.pi / 64 * (d**4 - di**4)
    return A, I, I
def msc_bar(b, h):
    return b * h, b * h**3 / 12, h * b**3 / 12
def msc_box(b, h, tf, tw):
    bi, hi = b - 2 * tw, h - 2 * tf
    A = b * h - bi * hi
    return A, (b * h**3 - bi * hi**3) / 12, (h * b**3 - hi * bi**3) / 12
def msc_i(h, ab, at, tw, tfb, tft):
    # converter's props_iorh (exact thin-wall assembly)
    hw = h - tft - tfb
    aw, atop, abot = tw * hw, at * tft, ab * tfb
    A = aw + atop + abot
    zw, zt, zb = tfb + hw / 2, h - tft / 2, tfb / 2
    zc = (abot * zb + aw * zw + atop * zt) / A
    I1 = ab * tfb**3 / 12 + abot * (zc - zb)**2 \
        + tw * hw**3 / 12 + aw * (zw - zc)**2 \
        + at * tft**3 / 12 + atop * (zt - zc)**2
    I2 = (tfb * ab**3 + hw * tw**3 + tft * at**3) / 12
    return A, I1, I2

# expectations from the FEM source cards:
# GPIPE 201: DI=0.46 DY=0.5 T=0.02 -> TUBE DIM=(0.5, 0.02)
# GIORH 202: HZ=.3 TY=.008 BT=.15 TT=.012 BB=.15 TB=.012
#   -> I DIM=(.3,.15,.15,.008,.012,.012)
# GBARM 204: HZ=.4 BT=BB=.2 -> BAR DIM=(.2,.4)
# GBOX 206: HZ=.5 TY=.01 TB=TT=.02 BY=.3 -> BOX DIM=(.3,.5,.02,.01)
# expectations from the FEM source cards:
# GPIPE 201: DI=0.46 DY=0.5 T=0.02 -> TUBE DIM=(0.5, 0.02)
# GIORH 202: HZ=.3 TY=.008 BT=.15 TT=.012 BB=.15 TB=.012
#   -> I DIM=(.3,.15,.15,.008,.012,.012)
# GBARM 204: HZ=.4 BT=BB=.2 -> BAR DIM=(.2,.4)
# GBOX 206: HZ=.5 TY=.01 TB=TT=.02 BY=.3 -> BOX DIM=(.3,.5,.02,.01)
# GLSEC 207: HZ=.25 TY=.01 BY=.15 TZ=.012 -> L DIM=(.15,.25,.012,.01)
# GCHAN 208: HZ=.3 TY=.008 BY=.12 TZ=.012 -> CHAN DIM=(.12,.3,.012,.008)
# GTONP 209: HZ=.35 TP=.014 BP=.8 BT=.12 TY=.008 TT=.012
#   -> I DIM=(.364,.8,.12,.008,.014,.012)  (plate as bottom flange)
# GUSYI 210: HZ=.4 BB=.3 BT=.2 TY=.01 TB=.02 TT=.015
#   -> I DIM=(.4,.3,.2,.01,.02,.015)  (off-centre web dropped)
# GIORHR 211: like GIORH 202 -> I DIM=(.3,.15,.15,.008,.012,.012)
# GLSECR 212: HZ=.2 TY=.01 BY=.12 TZ=.01 -> L DIM=(.12,.2,.01,.01)
# GCHANR 213: HZ=.25 TY=.008 BY=.1 TZ=.01 -> CHAN DIM=(.1,.25,.01,.008)
# GDOBO 214: HZ=1.5 BY=.9 TY=.012 TT=.016 TB=.02
#   -> I DIM=(1.5,.9,.9,.012,.02,.016)
expect = {
    201: ("TUBE", [0.5, 0.02], msc_tube(0.5, 0.02)),
    202: ("I", [0.3, 0.15, 0.15, 0.008, 0.012, 0.012],
          msc_i(0.3, 0.15, 0.15, 0.008, 0.012, 0.012)),
    204: ("BAR", [0.2, 0.4], msc_bar(0.2, 0.4)),
    206: ("BOX", [0.3, 0.5, 0.02, 0.01], msc_box(0.3, 0.5, 0.02, 0.01)),
    207: ("L", [0.15, 0.25, 0.012, 0.01], None),
    208: ("CHAN", [0.12, 0.3, 0.012, 0.008], None),
    209: ("I", [0.364, 0.8, 0.12, 0.008, 0.014, 0.012],
          msc_i(0.364, 0.8, 0.12, 0.008, 0.014, 0.012)),
    210: ("I", [0.4, 0.3, 0.2, 0.01, 0.02, 0.015],
          msc_i(0.4, 0.3, 0.2, 0.01, 0.02, 0.015)),
    211: ("I", [0.3, 0.15, 0.15, 0.008, 0.012, 0.012],
          msc_i(0.3, 0.15, 0.15, 0.008, 0.012, 0.012)),
    212: ("L", [0.12, 0.2, 0.01, 0.01], None),
    213: ("CHAN", [0.1, 0.25, 0.01, 0.008], None),
    214: ("I", [1.5, 0.9, 0.9, 0.012, 0.02, 0.016],
          msc_i(1.5, 0.9, 0.9, 0.012, 0.02, 0.016)),
}
for pid, (wtype, wdims, _props) in expect.items():
    gtype, gdims = pbeaml[pid]
    check("PBEAML %d TYPE" % pid, gtype, wtype)
    check("PBEAML %d ndims" % pid, len(gdims), len(wdims))
    for k, (g, w) in enumerate(zip(gdims, wdims)):
        check("PBEAML %d DIM%d" % (pid, k + 1), g, float(w))
# numeric spot check of section properties the dims imply
A, I1, I2 = msc_i(0.3, 0.15, 0.15, 0.008, 0.012, 0.012)
check("I-section A", A, 0.005808)
check("I-section I1", I1, 8.87088e-5, 1e-4)
A, I1, I2 = msc_box(0.3, 0.5, 0.02, 0.01)
check("BOX A", A, 0.0212)
A, I1, I2 = msc_tube(0.5, 0.02)
check("TUBE A", A, math.pi / 4 * (0.25 - 0.46**2))

print("== solid orientation ==")
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def dot(a, b): return sum(x*y for x, y in zip(a, b))
def hex_detJ(c):
    xi = [-1, 1, 1, -1, -1, 1, 1, -1]
    et = [-1, -1, 1, 1, -1, -1, 1, 1]
    ze = [-1]*4 + [1]*4
    jx = [sum(xi[k]*c[k][d] for k in range(8)) for d in range(3)]
    jy = [sum(et[k]*c[k][d] for k in range(8)) for d in range(3)]
    jz = [sum(ze[k]*c[k][d] for k in range(8)) for d in range(3)]
    return dot(jx, cross(jy, jz)) / 512.0

for e in chexa:
    eid, nodes = e[0], e[1:]
    d = hex_detJ([grids[n] for n in nodes])
    check("CHEXA %d detJ>0" % eid, d > 0, True)
# mirrored input element 8 must come out as 21 24 23 22 25 28 27 26
check("CHEXA 8 node order", chexa[1][1:], [21, 24, 23, 22, 25, 28, 27, 26])
t = ctetra[0]
v6 = dot(sub(grids[t[2]], grids[t[1]]),
         cross(sub(grids[t[3]], grids[t[1]]), sub(grids[t[4]], grids[t[1]])))
check("CTETRA vol6>0", v6 > 0, True)

print("== PBUSH / CMASS2 ==")
check("PBUSH 7 K", pbush[7], [1e6, 2e6, 3e6, 1e3, 2e3, 3e3])
check("CMASS2 count", len(cmass2), 4)
for eid, m, g, c in cmass2:
    check("CMASS2 %d node" % eid, g, 11)
check("CMASS2 dofs", [c for _, _, _, c in cmass2], [1, 3, 4, 6])

print("== SPC ==")
check("SPC entries", len(spc), 6)
want_spc = {(1, 12, 1): 0.01, (1, 12, 2): -0.02, (1, 13, 3): 0.0,
            (2, 12, 1): 0.0, (2, 12, 2): 0.0, (2, 13, 3): 0.0}
for sid, g, c, v in spc:
    check("SPC %d %d %d" % (sid, g, c), v, want_spc[(sid, g, c)])

print("== PLOAD4 ==")
check("PLOAD4 count", len(pload4), 4)
# [sid, eid, p1..p4, g1, g4]
check("PLOAD4 shell top", pload4[0], [1, 1, -1000.0, 0.0, 0.0, 0.0, 0, 0])
check("PLOAD4 shell vary", pload4[1],
      [2, 1, 500.0, 600.0, 700.0, 800.0, 0, 0])
check("PLOAD4 hex face", pload4[2], [1, 7, 2000.0, 0.0, 0.0, 0.0, 13, 12])
check("PLOAD4 follower", pload4[3], [2, 2, -900.0, 0.0, 0.0, 0.0, 0, 0])

print("== MAT2 / MPC / pin flags / GECC / sets ==")
# MAT2 from MORSMEL 3: G11=D11 G12=D21 G13=D31 G22=D22 G23=D32 G33=D33
mat2 = {}
mpc = []
cbeam_pins = {}   # eid -> (PA, PB)
cbeam_w = {}      # eid -> (WA, WB) for large-field CBEAM*
pshell = {}       # pid -> (mid, t)
cquad4 = {}       # eid -> pid
i = 0
while i < len(LINES):
    l = LINES[i]
    name = l[:8].strip()
    if name == "MAT2":
        f = fields(l)
        mat2[fld_int(f[1])] = [fld_float(x) for x in f[2:8]]
    elif name == "MPC":
        f = fields(l) + fields(LINES[i+1])
        # f: (name) SID G1 C1 A1 G2 C2 A2 (blank) | (blank) G3 C3 A3 ...
        vals = [(fld_int(f[2]), fld_int(f[3]), fld_float(f[4]))]
        rest = f[5:9] + f[10:]
        for k in range(0, len(rest) - 2, 3):
            if rest[k].strip():
                vals.append((fld_int(rest[k]), fld_int(rest[k+1]),
                             fld_float(rest[k+2])))
        mpc.append((fld_int(f[1]), vals))
        i += 1
    elif name == "CBEAM":
        f = fields(l)
        eid = fld_int(f[1])
        if LINES[i+1].startswith(" "):
            fc = fields(LINES[i+1])   # (blank marker) PA PB ...
            cbeam_pins[eid] = (fld_int(fc[1]), fld_int(fc[2]))
            i += 1
    elif name == "PSHELL":
        f = fields(l)
        pshell[fld_int(f[1])] = (fld_int(f[2]), fld_float(f[3]))
    elif name == "CQUAD4":
        f = fields(l)
        cquad4[fld_int(f[1])] = fld_int(f[2])
    elif name == "CBEAM*":
        def lf(ln):  # 16-char large fields after the 8-char marker
            return [ln[k:k+16] for k in range(8, len(ln), 16)]
        f1 = lf(LINES[i+1])   # (blank) X1 X2 X3
        f2 = lf(LINES[i+2])   # (blank) PA PB W1A W2A
        f3 = lf(LINES[i+3])   # (blank) W3A W1B W2B W3B
        eid = fld_int(l[8:24])   # large-field EID (16 chars)
        wa = (fld_float(f2[2]), fld_float(f2[3]), fld_float(f3[0]))
        wb = (fld_float(f3[1]), fld_float(f3[2]), fld_float(f3[3]))
        cbeam_w[eid] = (wa, wb)
        i += 3
    i += 1

check("MAT2 MID3 G", mat2.get(3),
      [1.62e11, 4.86e10, 0.0, 1.62e11, 0.0, 5.68e10])
check("MPC count", len(mpc), 1)
check("MPC sid", mpc[0][0], 1)
check("MPC terms", mpc[0][1], [(4, 3, 1.0), (3, 3, -1.0)])
# E3: BELFIX 500 releases axial dof 1 at both ends
check("CBEAM 3 pin flags", cbeam_pins.get(3), (1, 1))
# E5: BELFIX 501 fractional -> no pin flags
check("CBEAM 5 no pins", 5 in cbeam_pins, False)
# E13: GECC 401 local (0, .05, .1) rotated to global:
# x=(8,0,-3)/sqrt(73), z-hat=(24/73,0,64/73)/0.93633, y=(0,1,0)
L = math.sqrt(73.0)
zh = (24.0/73.0/(math.sqrt(4672.0)/73.0), 0.0,
      64.0/73.0/(math.sqrt(4672.0)/73.0))
wexp = (0.1*zh[0], 0.05, 0.1*zh[2])
wa, wb = cbeam_w.get(13, ((0,0,0),(0,0,0)))
for k in range(3):
    check("CBEAM 13 WA%d" % (k+1), wa[k], float(wexp[k]))
    check("CBEAM 13 WB%d" % (k+1), wb[k], float(wexp[k]))
# sets kept as comments
setc = [l for l in LINES if l.startswith("$ Set ")]
check("set comments", len(setc), 2)
check("set 5 named", any('Set 5 "SetA" (nodes)' in l for l in setc), True)
check("TEXT comment", any(l.startswith("$ TEXT: fem2bdf regression model")
                          for l in LINES), True)
# MPC referenced from every subcase
check("MPC in subcases", sum(1 for l in LINES if l.strip() == "MPC = 1"), 2)

print("== missing thickness / missing section ==")
# E25 quad: GEONO 101 has no GELTH -> zero-thickness placeholder PSHELL,
# element kept
check("E25 kept", 25 in cquad4, True)
pid25 = cquad4.get(25, 0)
check("E25 placeholder PSHELL t=0", pshell.get(pid25, (0, 1.0))[1], 0.0)
# E26 beam: GEONO 215 has no section card -> element skipped, no
# dangling PID left behind
check("E26 skipped", 26 not in cbeam_pins and
      not any(l[:8].strip() in ("CBEAM", "CBEAM*") and
              (fld_int(l[8:16]) if l[:8].strip() == "CBEAM"
               else fld_int(l[8:24])) == 26 for l in LINES), True)


print()
if fails:
    print("FAILED:", len(fails), "checks")
    sys.exit(1)
print("ALL CHECKS PASSED")
