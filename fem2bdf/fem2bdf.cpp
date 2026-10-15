/**
   \file fem2bdf.cpp
   \brief Convert a DNV GeniE / Sesam input interface file (FEM, T-file)
          into a Nastran Bulk Data File (BDF), built on top of DNV's
          open-source ExtFEMIO library (https://github.com/dnvgl/ExtFEMIO).

   Supported conversion (Sesam card -> BDF card):

     GNODE + GCOORD  -> GRID
     MISOSEL         -> MAT1
     GELTH           -> PSHELL
     GBEAMG          -> PBEAM   (A=AREA, I1=IY, I2=IZ, J=IX, I12=-IYZ,
                                 K1=SHARZ/AREA, K2=SHARY/AREA when non-zero)
     GPIPE           -> PBEAML (TYPE="TUBE", DIM=[DY, T])
     GIORH / GIORHR  -> PBEAML (TYPE="I", DIM=[HZ,BB,BT,TY,TB,TT];
                                 GIORHR fillet radii dropped)
     GUSYI           -> PBEAML (TYPE="I"; off-centre web dropped)
     GDOBO           -> PBEAML (TYPE="I", equal flanges BY, plates as flanges)
     GTONP           -> PBEAML (TYPE="I"; effective plate as bottom flange)
     GBOX            -> PBEAML (TYPE="BOX"/"BOX1")
     GCHAN / GCHANR  -> PBEAML (TYPE="CHAN", DIM=[BY,HZ,TZ,TY]; R dropped)
     GLSEC / GLSECR  -> PBEAML (TYPE="L",    DIM=[BY,HZ,TZ,TY]; R dropped)
     GBARM           -> PBEAML (TYPE="BAR"; trapezoidal section falls back
                                 to an equivalent PBEAM)
     FQUS/FFQ (24)   -> CQUAD4
     FTRS/FFTR (25)  -> CTRIA3
     BEAS (15) / BEPS (2) -> CBEAM (orientation vector = GUNIVEC, i.e. the
                                 Sesam local z-axis, which corresponds to
                                 the Nastran local y-axis; GECCEN -> WA/WB)
     TESS (10)       -> CROD
     LHEX (21) / IHEX (20)  -> CHEXA   (+ PSOLID)
     TETR (33) / ITET (31)  -> CTETRA  (+ PSOLID)
     TPRI (32) / IPRI (30)  -> CPENTA  (+ PSOLID)
     GSPR (18) + MGSPRNG    -> CBUSH to ground + PBUSH (diagonal terms)
     GMAS (11) + MGMASS     -> CONM1
     BNMASS          -> CMASS2 (one per dof)
     BNBCD (fixed)   -> SPC1 (written directly, ExtFEMIO has no SPC card)
     BNBCD (prescribed) + BNDISPL -> SPC with enforced values
     BNLOAD          -> FORCE / MOMENT (one load set per FEM load case)
     BEUSLO (normal pressure) -> PLOAD4 (shells and solid faces)
     TDLOAD          -> used for SUBCASE labels

   Everything else (dampers, follower/component-form surface loads,
   superelement hierarchy, ...) is reported and skipped.

   Usage:
     fem2bdf <input.fem> [-o <output.bdf>]
*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "fem/cards.h"
#include "fem/file.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "bdf/cards.h"
#include "bdf/cards_elements.h"
#include "bdf/cards_loads.h"
#include "bdf/cards_properties.h"

namespace fc = dnvgl::extfem::fem::cards;
namespace fe = dnvgl::extfem::fem::elements;
namespace bc = dnvgl::extfem::bdf::cards;

namespace {

/// Collected content of the FEM input file.
struct FemModel {
    // nodes: internal number -> external number / coordinates / dof list
    std::map<long, long> node_ext;              ///< NODENO -> NODEX
    std::map<long, fc::gcoord> coords;          ///< NODENO -> GCOORD
    std::map<long, std::vector<int>> node_dofs; ///< NODENO -> ODOF

    std::map<long, fc::gelmnt1> elems;          ///< ELNO -> GELMNT1
    std::map<long, fc::gelref1> refs;           ///< ELNO -> GELREF1

    std::map<long, fc::misosel> mats;           ///< MATNO -> MISOSEL

    std::map<long, fc::gelth>   thickness;      ///< GEONO -> GELTH
    std::map<long, fc::gbeamg>  beamg;          ///< GEONO -> GBEAMG
    std::map<long, fc::gpipe>   pipes;          ///< GEONO -> GPIPE
    std::map<long, fc::giorh>   iorh;           ///< GEONO -> GIORH
    std::map<long, fc::giorhr>  iorhrs;         ///< GEONO -> GIORHR
    std::map<long, fc::gusyi>   usyis;          ///< GEONO -> GUSYI
    std::map<long, fc::gbarm>   barm;           ///< GEONO -> GBARM
    std::map<long, fc::gbox>    boxes;          ///< GEONO -> GBOX
    std::map<long, fc::gdobo>   dobos;          ///< GEONO -> GDOBO
    std::map<long, fc::gchan>   chans;          ///< GEONO -> GCHAN
    std::map<long, fc::gchanr>  chanrs;         ///< GEONO -> GCHANR
    std::map<long, fc::glsec>   lsecs;          ///< GEONO -> GLSEC
    std::map<long, fc::glsecr>  lsecrs;         ///< GEONO -> GLSECR
    std::map<long, fc::gtonp>   tonps;          ///< GEONO -> GTONP

    std::map<long, fc::geccen>  eccs;           ///< ECCNO -> GECCEN (global)
    std::map<long, fc::gecc>    leccs;          ///< ECCNO -> GECC (element local)
    std::map<long, fc::gunivec> trans;          ///< TRANSNO -> GUNIVEC
    std::map<long, fc::belfix>  hinges;         ///< FIXNO -> BELFIX
    std::vector<fc::bldep>  ldeps;              ///< BLDEP linear dependences

    std::map<long, fc::mgsprng> springs;        ///< MATNO -> MGSPRNG
    std::map<long, fc::morsmel> morsmels;       ///< MATNO -> MORSMEL
    std::map<long, std::vector<double>> mgmass; ///< MATNO -> MGMASS (lower tri)
    struct NodalMass {
        long nodeno;
        std::vector<double> mass;
    };
    std::vector<NodalMass> nodal_masses;        ///< BNMASS cards

    /// Sesam sets: ISREF -> (type 1=node/2=elem, members accumulated
    /// across GSETMEMB records); names come from TDSETNAM.
    std::map<long, std::pair<long, std::vector<long>>> sets;
    std::map<long, std::string> set_names;      ///< ISREF -> TDSETNAM name
    std::vector<std::string> user_texts;        ///< TEXT card contents

    std::vector<fc::bnbcd>  bcs;                ///< BNBCD cards
    std::vector<fc::bndispl> disps;             ///< BNDISPL cards
    std::vector<fc::bnload> loads;              ///< BNLOAD cards
    std::vector<fc::beuslo> sloads;             ///< BEUSLO surface loads
    std::map<long, std::string> load_names;     ///< LLC -> TDLOAD name
};

/// Counters for the final report.
struct Stats {
    size_t cards_read{0};
    size_t cards_failed{0};
    std::map<std::string, size_t> skipped_cards;
    std::map<long, size_t> skipped_eltypes;
    size_t missing_geono{0};
    size_t no_geono{0};
    size_t missing_mat{0};
    size_t missing_spring{0};      ///< GSPR without MGSPRNG
    size_t missing_mass{0};        ///< GMAS without MGMASS
    size_t spring_coupling{0};     ///< springs with dropped off-diagonal K
    size_t spring_local{0};        ///< springs/masses w/ ignored local system
    size_t prescribed_bc{0};
    size_t missing_prescribed{0};  ///< FIX=2 but no BNDISPL for a case
    size_t accel_disp{0};          ///< BNDISPL DTYPE=3 skipped
    size_t complx_disp{0};
    size_t missing_elem_load{0};   ///< BEUSLO referencing unknown element
    size_t comp_form_load{0};      ///< BEUSLO LOTYP=2 skipped
    size_t hydro_dummy{0};         ///< BEUSLO LOTYP=3 skipped
    size_t follower_load{0};       ///< BEUSLO LOTYP=-1 (converted as normal)
    size_t complx_load{0};
    size_t side_mid{0};            ///< BEUSLO SIDE=2 (treated as SIDE=3)
    size_t bad_side{0};            ///< BEUSLO with unusable SIDE flag
    size_t trapezoid_bar{0};       ///< GBARM trapezoid -> equivalent PBEAM
    size_t usyi_offset{0};         ///< GUSYI with off-centre web (dropped)
    size_t curvature_ignored{0};   ///< GIORHR/GLSECR/GCHANR (R dropped)
    size_t tonp_as_i{0};           ///< GTONP converted to I section
    size_t default_mats{0};        ///< MATNO without MISOSEL; default MAT1 emitted
    size_t shcen_adjust{0};        ///< beam offsets shifted centroid->shear centre
    size_t gecc_local{0};          ///< GECC local eccentricities rotated to global
    size_t hinge_applied{0};       ///< beam ends with pin flags from BELFIX
    size_t hinge_partial{0};       ///< fractional fixation (0<A<1) treated as fixed
    size_t hinge_spring{0};        ///< OPT=SPRING hinge stiffness dropped
    size_t hinge_local{0};         ///< BELFIX with TRANO ignored (global assumed)
    size_t morsmel_q{0};           ///< MORSMEL direction vector Q dropped
    size_t missing_ldep{0};        ///< dependent dof (BNBCD=3) without BLDEP
    size_t missing_hinge{0};       ///< element references unknown BELFIX FIXNO
    size_t n_mpc{0};               ///< MPC cards written from BLDEP
    size_t default_thk{0};         ///< shells w/o GELTH (or TH<=0); default t emitted
    size_t no_prop{0};             ///< elements skipped: property could not be built
};

void warn_skipped(Stats &st, const std::string &name) {
    st.skipped_cards[name]++;
}

/// Parse an integer entry that may be written in float format
/// ("1.10000000E+01"); std::stol would stop at the decimal point.
long field_to_long(const std::string &s) {
    return std::lround(std::stod(s));
}

/// Manual fallback parser for BNLOAD cards written in the *documented*
/// layout (no padding field after COMPLX, i.e. NODENO at index 4).
/// ExtFEMIO's bnload reader expects the de-facto GeniE layout where the
/// first record line is padded and NODENO sits at index 5; real GeniE
/// files follow that convention, but hand-written files may not.
bool parse_bnload_official(const std::vector<std::string> &entries,
                           size_t nfields, FemModel &model) {
    try {
        if (nfields < 6) return false;
        const long llc = field_to_long(entries.at(1));
        const long lotyp = field_to_long(entries.at(2));
        const bool complx = std::stod(entries.at(3)) != 0.0;
        const long nodeno = field_to_long(entries.at(4));
        const long ndof = field_to_long(entries.at(5));
        if (ndof <= 0 || ndof > 12) return false;
        std::vector<double> rload, iload;
        for (long i = 0; i < ndof; ++i)
            rload.push_back(std::stod(entries.at(6 + i)));
        if (complx)
            for (long i = 0; i < ndof; ++i)
                iload.push_back(std::stod(entries.at(6 + ndof + i)));
        model.loads.emplace_back(
            fc::bnload(llc, lotyp, complx, nodeno, rload, iload));
        return true;
    } catch (const std::exception &) {
        return false;
    }
}

/// Manual parser for BNMASS (no ExtFEMIO reader):
///   BNMASS NODENO NDOF MASS1 ... MASS(NDOF)
bool parse_bnmass(const std::vector<std::string> &entries,
                  size_t nfields, FemModel &model) {
    try {
        if (nfields < 4) return false;
        FemModel::NodalMass nm;
        nm.nodeno = field_to_long(entries.at(1));
        const long ndof = field_to_long(entries.at(2));
        if (ndof <= 0 || ndof > 12 || nfields < 3 + (size_t)ndof) return false;
        for (long i = 0; i < ndof; ++i)
            nm.mass.push_back(std::stod(entries.at(3 + i)));
        model.nodal_masses.push_back(nm);
        return true;
    } catch (const std::exception &) {
        return false;
    }
}

/// Manual parser for MGMASS (no ExtFEMIO reader):
///   MGMASS MATNO NDOF M(1,1) M(2,1) ... M(n,1) M(2,2) ...
///   (lower triangle, column by column)
bool parse_mgmass(const std::vector<std::string> &entries,
                  size_t nfields, FemModel &model) {
    try {
        if (nfields < 4) return false;
        const long matno = field_to_long(entries.at(1));
        const long ndof = field_to_long(entries.at(2));
        if (ndof <= 0 || ndof > 12) return false;
        const size_t nvals = (size_t)(ndof * (ndof + 1) / 2);
        if (nfields < 3 + nvals) return false;
        std::vector<double> m;
        for (size_t i = 0; i < nvals; ++i)
            m.push_back(std::stod(entries.at(3 + i)));
        model.mgmass[matno] = m;
        return true;
    } catch (const std::exception &) {
        return false;
    }
}

/// Read the whole FEM file into a FemModel.
bool read_fem(const std::string &path, FemModel &model, Stats &st) {
    std::ifstream fs(path.c_str());
    if (!fs) {
        std::cerr << "ERROR: cannot open input file '" << path << "'\n";
        return false;
    }

    dnvgl::extfem::fem::input::fem_file inp(fs);
    std::vector<std::string> lines;
    std::vector<std::string> entries;
    std::unique_ptr<fc::__base::card> cur;

    while (!inp.eof()) {
        size_t len = inp.get(lines);
        if (len == 0) break;
        size_t nfields = fc::__base::card::card_split(lines, len, entries);
        st.cards_read++;

        // Cards without an ExtFEMIO reader are parsed manually.
        {
            const std::string cname = entries.empty() ? "" : entries.front();
            if (cname == "BNMASS") {
                if (!parse_bnmass(entries, nfields, model))
                    warn_skipped(st, "BNMASS (parse error)");
                continue;
            }
            if (cname == "MGMASS") {
                if (!parse_mgmass(entries, nfields, model))
                    warn_skipped(st, "MGMASS (parse error)");
                continue;
            }
        }

        try {
            fc::dispatch(entries, nfields, cur);
        } catch (const std::exception &) {
            std::string name = entries.empty() ? "?" : entries.front();
            // retry BNLOAD with the documented (unpadded) field layout
            if (name == "BNLOAD" &&
                parse_bnload_official(entries, nfields, model))
                continue;
            st.cards_failed++;
            warn_skipped(st, name + " (parse error)");
            continue;
        }

        try {
            switch (cur->card_type()) {
            case fc::types::GNODE: {
                auto *c = static_cast<fc::gnode*>(cur.get());
                model.node_ext[c->NODENO] = c->NODEX;
                model.node_dofs[c->NODENO] = c->ODOF;
                break;
            }
            case fc::types::GCOORD: {
                auto *c = static_cast<fc::gcoord*>(cur.get());
                model.coords[c->NODENO] = *c;
                break;
            }
            case fc::types::GELMNT1: {
                auto *c = static_cast<fc::gelmnt1*>(cur.get());
                model.elems[c->ELNO] = *c;
                break;
            }
            case fc::types::GELREF1: {
                auto *c = static_cast<fc::gelref1*>(cur.get());
                model.refs[c->ELNO] = *c;
                break;
            }
            case fc::types::MISOSEL: {
                auto *c = static_cast<fc::misosel*>(cur.get());
                model.mats[c->MATNO] = *c;
                break;
            }
            case fc::types::GELTH: {
                auto *c = static_cast<fc::gelth*>(cur.get());
                model.thickness[c->GEONO] = *c;
                break;
            }
            case fc::types::GBEAMG: {
                auto *c = static_cast<fc::gbeamg*>(cur.get());
                model.beamg[c->GEONO] = *c;
                break;
            }
            case fc::types::GPIPE: {
                auto *c = static_cast<fc::gpipe*>(cur.get());
                model.pipes[c->GEONO] = *c;
                break;
            }
            case fc::types::GIORH: {
                auto *c = static_cast<fc::giorh*>(cur.get());
                model.iorh[c->GEONO] = *c;
                break;
            }
            case fc::types::GIORHR: {
                auto *c = static_cast<fc::giorhr*>(cur.get());
                model.iorhrs[c->GEONO] = *c;
                break;
            }
            case fc::types::GUSYI: {
                auto *c = static_cast<fc::gusyi*>(cur.get());
                model.usyis[c->GEONO] = *c;
                break;
            }
            case fc::types::GBARM: {
                auto *c = static_cast<fc::gbarm*>(cur.get());
                model.barm[c->GEONO] = *c;
                break;
            }
            case fc::types::GBOX: {
                auto *c = static_cast<fc::gbox*>(cur.get());
                model.boxes[c->GEONO] = *c;
                break;
            }
            case fc::types::GDOBO: {
                auto *c = static_cast<fc::gdobo*>(cur.get());
                model.dobos[c->GEONO] = *c;
                break;
            }
            case fc::types::GCHAN: {
                auto *c = static_cast<fc::gchan*>(cur.get());
                model.chans[c->GEONO] = *c;
                break;
            }
            case fc::types::GCHANR: {
                auto *c = static_cast<fc::gchanr*>(cur.get());
                model.chanrs[c->GEONO] = *c;
                break;
            }
            case fc::types::GLSEC: {
                auto *c = static_cast<fc::glsec*>(cur.get());
                model.lsecs[c->GEONO] = *c;
                break;
            }
            case fc::types::GLSECR: {
                auto *c = static_cast<fc::glsecr*>(cur.get());
                model.lsecrs[c->GEONO] = *c;
                break;
            }
            case fc::types::GTONP: {
                auto *c = static_cast<fc::gtonp*>(cur.get());
                model.tonps[c->GEONO] = *c;
                break;
            }
            case fc::types::GECCEN: {
                auto *c = static_cast<fc::geccen*>(cur.get());
                model.eccs[c->ECCNO] = *c;
                break;
            }
            case fc::types::GECC: {
                // local-frame eccentricity (ExtFEMIO's doxygen wrongly
                // says superelement system; Sesam docs say local)
                auto *c = static_cast<fc::gecc*>(cur.get());
                model.leccs[c->ECCNO] = *c;
                break;
            }
            case fc::types::BELFIX: {
                auto *c = static_cast<fc::belfix*>(cur.get());
                model.hinges[c->FIXNO] = *c;
                break;
            }
            case fc::types::BLDEP: {
                model.ldeps.push_back(*static_cast<fc::bldep*>(cur.get()));
                break;
            }
            case fc::types::MORSMEL: {
                auto *c = static_cast<fc::morsmel*>(cur.get());
                model.morsmels[c->MATNO] = *c;
                break;
            }
            case fc::types::GSETMEMB: {
                auto *c = static_cast<fc::gsetmemb*>(cur.get());
                auto &s = model.sets[c->ISREF];
                s.first = (c->ISTYPE == fc::gsetmemb::types::ELEM_SET) ? 2 : 1;
                s.second.insert(s.second.end(),
                                c->IRMEMB.begin(), c->IRMEMB.end());
                break;
            }
            case fc::types::TDSETNAM: {
                auto *c = static_cast<fc::tdsetnam*>(cur.get());
                if (!c->SET_NAME.empty())
                    model.set_names[c->ISREF] = c->SET_NAME;
                break;
            }
            case fc::types::TEXT: {
                auto *c = static_cast<fc::text*>(cur.get());
                model.user_texts.insert(model.user_texts.end(),
                                        c->CONT.begin(), c->CONT.end());
                break;
            }
            case fc::types::GUNIVEC: {
                auto *c = static_cast<fc::gunivec*>(cur.get());
                model.trans[c->TRANSNO] = *c;
                break;
            }
            case fc::types::MGSPRNG: {
                auto *c = static_cast<fc::mgsprng*>(cur.get());
                model.springs[c->MATNO] = *c;
                break;
            }
            case fc::types::BNBCD: {
                model.bcs.push_back(*static_cast<fc::bnbcd*>(cur.get()));
                break;
            }
            case fc::types::BNDISPL: {
                model.disps.push_back(*static_cast<fc::bndispl*>(cur.get()));
                break;
            }
            case fc::types::BEUSLO: {
                model.sloads.push_back(*static_cast<fc::beuslo*>(cur.get()));
                break;
            }
            case fc::types::BNLOAD: {
                auto *c = static_cast<fc::bnload*>(cur.get());
                // Guard against silently mis-padded input: if the card
                // claims zero dofs, retry with the documented layout.
                if (c->NDOF == 0 &&
                    parse_bnload_official(entries, nfields, model))
                    break;
                model.loads.push_back(*c);
                break;
            }
            case fc::types::TDLOAD: {
                auto *c = static_cast<fc::tdload*>(cur.get());
                if (!c->SET_NAME.empty())
                    model.load_names[c->ILREF] = c->SET_NAME;
                break;
            }
            default:
                // DATE, IDENT, IEND, HIERARCH, HSUPSTAT, HSUPTRAN,
                // BSELL, GELMNT2, UNKNOWN, ...
                break;
            }
        } catch (const std::exception &) {
            st.cards_failed++;
        }
    }
    return true;
}

/// External node number for an internal one (1:1 fallback).
long ext_node(const FemModel &m, long nodeno) {
    auto it = m.node_ext.find(nodeno);
    return it != m.node_ext.end() ? it->second : nodeno;
}

/// Dof list of a node (defaults to 1..6).
std::vector<int> dofs_of(const FemModel &m, long nodeno, long ndof) {
    auto it = m.node_dofs.find(nodeno);
    if (it != m.node_dofs.end() && !it->second.empty()) return it->second;
    std::vector<int> d;
    for (int i = 1; i <= (ndof > 0 ? ndof : 6); ++i) d.push_back(i);
    return d;
}

/// GEONO referenced by an element through its GELREF1 (0 = none).
long elem_geono(const fc::gelref1 &r) {
    if (r.GEONO_OPT > 0) return r.GEONO_OPT;
    if (r.GEONO_OPT == -1 && !r.GEONO.empty()) return r.GEONO.front();
    return 0;
}

/// TRANSNO referenced by an element (first node), 0 = none.
long elem_transno(const fc::gelref1 &r) {
    if (r.TRANSNO_OPT > 0) return r.TRANSNO_OPT;
    if (r.TRANSNO_OPT == -1 && !r.TRANSNO.empty()) return r.TRANSNO.front();
    return 0;
}

/// True when GEONO resolves to a shape card written as PBEAML.
/// GBARM is excluded: rectangular bars have SHCEN=0 anyway, and
/// trapezoidal bars fall back to PBEAM (centroidal axis).
bool shape_card_for(const FemModel &m, long geono) {
    return m.pipes.count(geono) || m.iorh.count(geono) ||
           m.iorhrs.count(geono) || m.usyis.count(geono) ||
           m.dobos.count(geono) || m.tonps.count(geono) ||
           m.chans.count(geono) || m.chanrs.count(geono) ||
           m.lsecs.count(geono) || m.lsecrs.count(geono) ||
           m.boxes.count(geono);
}

/// Computed section properties for cross-section cards (used as
/// fallback when a shape cannot be mapped to a native PBEAML type).
struct SectionProps {
    double area{0}, iy{0}, iz{0}, j{0};
};

/// GIORH (I or H section): HZ total height, TY web thickness,
/// BT/TT top flange width/thickness, BB/TB bottom flange.
SectionProps props_iorh(const fc::giorh &s) {
    const double hw = std::max(0.0, s.HZ - s.TT - s.TB);
    const double a_web = s.TY * hw;
    const double a_top = s.BT * s.TT;
    const double a_bot = s.BB * s.TB;

    SectionProps p;
    p.area = a_web + a_top + a_bot;
    if (p.area <= 0.0) return p;

    // centroid measured from the bottom fibre
    const double z_web = s.TB + hw / 2.0;
    const double z_top = s.HZ - s.TT / 2.0;
    const double z_bot = s.TB / 2.0;
    const double zc =
        (a_bot * z_bot + a_web * z_web + a_top * z_top) / p.area;

    // IY = integral z^2 dA (strong axis for an upright I)
    p.iy = s.BB * std::pow(s.TB, 3) / 12.0 + a_bot * std::pow(zc - z_bot, 2)
         + s.TY * std::pow(hw, 3) / 12.0 + a_web * std::pow(z_web - zc, 2)
         + s.BT * std::pow(s.TT, 3) / 12.0 + a_top * std::pow(z_top - zc, 2);
    // IZ = integral y^2 dA
    p.iz = (s.TB * std::pow(s.BB, 3) + hw * std::pow(s.TY, 3)
            + s.TT * std::pow(s.BT, 3)) / 12.0;
    // St. Venant torsion constant (thin-walled approximation)
    p.j = (s.BB * std::pow(s.TB, 3) + hw * std::pow(s.TY, 3)
           + s.BT * std::pow(s.TT, 3)) / 3.0;
    return p;
}

/// Torsion constant of a solid rectangle (b >= t).
double torsion_rect(double b, double t) {
    if (b < t) std::swap(b, t);
    if (b <= 0.0 || t <= 0.0) return 0.0;
    const double r = t / b;
    return b * std::pow(t, 3) * (1.0 / 3.0 - 0.21 * r * (1.0 - std::pow(r, 4) / 12.0));
}

/// GBARM (massive bar, possibly trapezoidal): HZ height,
/// BT top width, BB bottom width.
SectionProps props_gbarm(const fc::gbarm &s) {
    SectionProps p;
    const double a = s.BB, b = s.BT, h = s.HZ;
    p.area = (a + b) * h / 2.0;
    if (p.area <= 0.0) return p;
    // exact trapezoid inertia about the centroidal y axis (bases horizontal)
    p.iy = std::pow(h, 3) * (a * a + 4.0 * a * b + b * b) / (36.0 * (a + b));
    // exact symmetric-trapezoid inertia about the vertical centroidal axis
    p.iz = h * (std::pow(a, 3) + a * a * b + a * b * b + std::pow(b, 3)) / 48.0;
    // torsion: rectangle approximation using the mean width
    p.j = torsion_rect(h, (a + b) / 2.0);
    return p;
}

/// GBOX (box girder): BY width, HZ height, TY web thickness (2 webs),
/// TT top / TB bottom flange thickness.
SectionProps props_gbox(const fc::gbox &s) {
    SectionProps p;
    const double h = s.HZ, b = s.BY, tw = s.TY;
    const double tt = s.TT, tb = s.TB;
    const double hi = h - tt - tb;          // clear inner height
    const double bi = b - 2.0 * tw;         // clear inner width
    p.area = b * h - bi * hi;
    if (p.area <= 0.0) return p;
    p.iy = (b * std::pow(h, 3) - bi * std::pow(hi, 3)) / 12.0;
    p.iz = (h * std::pow(b, 3) - hi * std::pow(bi, 3)) / 12.0;
    // Bredt-Batho torsion of the closed cell (median line)
    const double lv = h - (tt + tb) / 2.0;  // median height
    const double lh = b - tw;               // median width
    const double am = lv * lh;              // enclosed median area
    double ds_t = 0.0;
    if (tw > 0) ds_t += 2.0 * lv / tw;
    if (tt > 0) ds_t += lh / tt;
    if (tb > 0) ds_t += lh / tb;
    p.j = ds_t > 0.0 ? 4.0 * am * am / ds_t : 0.0;
    return p;
}

/// Pipe area from GPIPE data.
double pipe_area(const fc::gpipe &s) {
    const double t = s.T > 0.0 ? s.T : (s.DY - s.DI) / 2.0;
    return M_PI / 4.0 * (s.DY * s.DY - std::pow(s.DY - 2.0 * t, 2));
}

/// PID allocator: keeps property ids unique across all property types
/// (ExtFEMIO enforces global uniqueness of PID).
class PidPool {
public:
    explicit PidPool(long start) : next_free_(start) {}

    /// Get (or create) the PID for a (GEONO, MATNO) pair. The first
    /// material seen for a GEONO keeps the GEONO value as PID.
    long get(long geono, long matno) {
        const auto key = std::make_pair(geono, matno);
        auto it = map_.find(key);
        if (it != map_.end()) return it->second;
        long pid;
        if (seen_geono_.insert(geono).second && used_.insert(geono).second) {
            pid = geono;
        } else {
            while (used_.count(next_free_) || map_.count(std::make_pair(next_free_, next_free_)))
                ++next_free_;
            pid = next_free_;
            used_.insert(pid);
            ++next_free_;
        }
        map_[key] = pid;
        return pid;
    }

private:
    long next_free_;
    std::set<long> used_;
    std::set<long> seen_geono_;
    std::map<std::pair<long, long>, long> map_;
};

const char *eltype_name(fe::el_types t) {
    switch (t) {
    case fe::el_types::BEPS: return "BEPS";
    case fe::el_types::TESS: return "TESS";
    case fe::el_types::BEAS: return "BEAS";
    case fe::el_types::FQUS_FFQ: return "FQUS/FFQ";
    case fe::el_types::FTRS_FFTR: return "FTRS/FFTR";
    case fe::el_types::BTSS: return "BTSS";
    case fe::el_types::GMAS: return "GMAS";
    case fe::el_types::GSPR: return "GSPR";
    case fe::el_types::GLSH: return "GLSH";
    case fe::el_types::LHEX: return "LHEX";
    case fe::el_types::IHEX: return "IHEX";
    case fe::el_types::TETR: return "TETR";
    case fe::el_types::ITET: return "ITET";
    case fe::el_types::TPRI: return "TPRI";
    case fe::el_types::IPRI: return "IPRI";
    case fe::el_types::LQUA: return "LQUA";
    case fe::el_types::ILST: return "ILST";
    case fe::el_types::IQQE: return "IQQE";
    default: return "other";
    }
}

bool is_solid(fe::el_types t) {
    switch (t) {
    case fe::el_types::LHEX:
    case fe::el_types::IHEX:
    case fe::el_types::TETR:
    case fe::el_types::ITET:
    case fe::el_types::TPRI:
    case fe::el_types::IPRI:
        return true;
    default:
        return false;
    }
}

// ------------------------------------------------------------------
// Minimal fixed-format helpers for the few BDF cards ExtFEMIO does
// not implement (PSOLID, CHEXA/CTETRA/CPENTA, SPC, PLOAD4).
// ------------------------------------------------------------------

/// Integer as an 8-char field.
std::string f8i(long v) {
    char buf[24];
    std::snprintf(buf, sizeof(buf), "%8ld", v);
    return buf;
}

/// Real as an 8-char field. Uses plain decimal with as many digits as
/// fit; falls back to the Nastran exponent form without 'E' ("1.00+06")
/// which always fits 8 chars. (printf %8.4f would overflow the field
/// for |v| >= 1000 and corrupt all following fields.)
std::string f8r(double v) {
    char buf[32];
    if (v == 0.0) {
        std::snprintf(buf, sizeof(buf), "%8.4f", 0.0);
        return buf;
    }
    const double a = std::fabs(v);
    if (a >= 1.0e-3 && a < 1.0e6) {
        for (int prec = 4; prec >= 0; --prec) {
            std::snprintf(buf, sizeof(buf), "%8.*f", prec, v);
            if (std::strlen(buf) <= 8) return buf;
        }
    }
    int e = (int)std::floor(std::log10(a));
    double mant = v / std::pow(10.0, (double)e);
    if (std::fabs(mant) >= 9.995) { mant /= 10.0; ++e; }
    for (int prec = 2; prec >= 1; --prec) {
        std::snprintf(buf, sizeof(buf), "%.*f%c%02d", prec, mant,
                      e < 0 ? '-' : '+', std::abs(e));
        if (std::strlen(buf) <= 8) {
            const std::string r = buf;
            return std::string(8 - r.size(), ' ') + r;
        }
    }
    std::snprintf(buf, sizeof(buf), "%8.0E", v);  // last resort
    return buf;
}

// ------------------------------------------------------------------
// Solid element geometry: node-order orientation check / fix.
// ------------------------------------------------------------------

struct Vec3 {
    double x{0}, y{0}, z{0};
};

Vec3 sub(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}
double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

/// Coordinates of the element corner nodes (internal numbering).
std::vector<Vec3> corner_coords(const FemModel &m,
                                const std::vector<long> &nodin,
                                size_t ncorner) {
    std::vector<Vec3> c;
    for (size_t i = 0; i < ncorner && i < nodin.size(); ++i) {
        auto it = m.coords.find(nodin[i]);
        if (it == m.coords.end()) return {};
        c.push_back({it->second.XCOORD, it->second.YCOORD,
                     it->second.ZCOORD});
    }
    return c;
}

/// det J at the centre of an 8-node brick (iso corner order:
/// 1=(-1,-1,-1) 2=(+1,-1,-1) 3=(+1,+1,-1) 4=(-1,+1,-1), 5..8 top).
double hex_detJ(const std::vector<Vec3> &c) {
    static const int xi[8]  = {-1, 1, 1, -1, -1, 1, 1, -1};
    static const int eta[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
    static const int zet[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    Vec3 jx, jy, jz;
    for (int i = 0; i < 8; ++i) {
        jx.x += xi[i] * c[i].x; jx.y += xi[i] * c[i].y; jx.z += xi[i] * c[i].z;
        jy.x += eta[i] * c[i].x; jy.y += eta[i] * c[i].y; jy.z += eta[i] * c[i].z;
        jz.x += zet[i] * c[i].x; jz.y += zet[i] * c[i].y; jz.z += zet[i] * c[i].z;
    }
    return dot(jx, cross(jy, jz)) / 512.0;
}

/// Signed volume x6 of a tetrahedron.
double tet_vol6(const std::vector<Vec3> &c) {
    return dot(sub(c[1], c[0]), cross(sub(c[2], c[0]), sub(c[3], c[0])));
}

/// det J at the centre of a 6-node wedge
/// (1=(0,0,-1) 2=(1,0,-1) 3=(0,1,-1), 4..6 top at t=+1).
double wedge_detJ(const std::vector<Vec3> &c) {
    // J_r = 0.5*(x2-x1) + 0.5*(x5-x4);  J_s analogous; J_t = (top-bottom)/6
    Vec3 jr{0.5 * (c[1].x - c[0].x) + 0.5 * (c[4].x - c[3].x),
            0.5 * (c[1].y - c[0].y) + 0.5 * (c[4].y - c[3].y),
            0.5 * (c[1].z - c[0].z) + 0.5 * (c[4].z - c[3].z)};
    Vec3 js{0.5 * (c[2].x - c[0].x) + 0.5 * (c[5].x - c[3].x),
            0.5 * (c[2].y - c[0].y) + 0.5 * (c[5].y - c[3].y),
            0.5 * (c[2].z - c[0].z) + 0.5 * (c[5].z - c[3].z)};
    Vec3 jt;
    for (int i = 0; i < 3; ++i) {
        jt.x += (c[3 + i].x - c[i].x) / 6.0;
        jt.y += (c[3 + i].y - c[i].y) / 6.0;
        jt.z += (c[3 + i].z - c[i].z) / 6.0;
    }
    return dot(jr, cross(js, jt));
}

// ------------------------------------------------------------------
// Solid face tables for BEUSLO -> PLOAD4.
// Sesam local corner numbers (1-based) per side, ordered so the face
// normal points outwards (checked geometrically at run time).
// ------------------------------------------------------------------

/// LHEX / IHEX: 6 quad faces; g4[s] = corner adjacent to face-node 1
/// that is NOT on the face (PLOAD4 face identification).
const int HEX_FACE[6][4] = {
    {1, 2, 3, 4}, {1, 2, 6, 5}, {2, 3, 7, 6},
    {3, 4, 8, 7}, {4, 1, 5, 8}, {5, 6, 7, 8}
};
const int HEX_G4[6] = {5, 5, 1, 2, 3, 1};

/// TETR / ITET: 4 triangular faces; g4[s] = opposite corner.
const int TET_FACE[4][3] = {
    {3, 2, 4}, {1, 3, 4}, {2, 1, 4}, {1, 2, 3}
};
const int TET_G4[4] = {1, 2, 3, 4};

/// TPRI / IPRI: sides 1-3 quad, 4-5 triangular; g4 for the quad faces.
const int PRI_FACE[5][4] = {
    {2, 5, 3, 6}, {3, 6, 1, 4}, {1, 4, 2, 5},
    {1, 2, 3, 0}, {6, 5, 4, 0}
};
const int PRI_G4[5] = {1, 2, 3, 0, 0};

void usage(const char *prog) {
    std::cerr
        << "fem2bdf - convert DNV GeniE/Sesam FEM (T-file) to Nastran BDF\n"
        << "          built on ExtFEMIO (https://github.com/dnvgl/ExtFEMIO)\n\n"
        << "usage: " << prog << " <input.fem> [-o <output.bdf>]\n";
}

} // anonymous namespace

int main(int argc, char **argv) {
    std::string in_path, out_path;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-o" && i + 1 < argc) {
            out_path = argv[++i];
        } else if (a == "-h" || a == "--help") {
            usage(argv[0]);
            return 0;
        } else if (in_path.empty()) {
            in_path = a;
        } else {
            usage(argv[0]);
            return 2;
        }
    }
    if (in_path.empty()) {
        usage(argv[0]);
        return 2;
    }
    if (out_path.empty()) {
        out_path = in_path;
        const auto dot = out_path.find_last_of('.');
        const auto sep = out_path.find_last_of("/\\");
        if (dot != std::string::npos &&
            (sep == std::string::npos || dot > sep))
            out_path.erase(dot);
        out_path += ".bdf";
    }

    FemModel model;
    Stats st;
    if (!read_fem(in_path, model, st)) return 1;

    std::ofstream out(out_path.c_str());
    if (!out) {
        std::cerr << "ERROR: cannot open output file '" << out_path << "'\n";
        return 1;
    }

    // --------------------------------------------------------------
    // Load cases: nodal loads, surface loads and prescribed
    // displacements all carry a load-case number.
    // --------------------------------------------------------------
    std::set<long> load_cases;
    for (const auto &l : model.loads)  load_cases.insert(l.LLC);
    for (const auto &l : model.sloads) load_cases.insert(l.LLC);
    for (const auto &d : model.disps)  load_cases.insert(d.LLC);

    // Map constrained-dof pattern -> external node numbers (for SPC1)
    std::map<std::string, std::vector<long>> spc_groups;
    bool any_prescribed_flag = false;
    for (const auto &b : model.bcs) {
        const auto dofs = dofs_of(model, b.NODENO, b.NDOF);
        std::string comp;
        bool prescribed = false;
        for (size_t i = 0; i < b.FIX.size() && i < dofs.size(); ++i) {
            switch (b.FIX[i]) {
            case fc::bnbcd::fix_key::DISPL_FIX:
                comp += char('0' + dofs[i]);
                break;
            case fc::bnbcd::fix_key::PRESCRIBED:
                prescribed = true;
                break;
            default:
                break;
            }
        }
        if (prescribed) { st.prescribed_bc++; any_prescribed_flag = true; }
        if (!comp.empty())
            spc_groups[comp].push_back(ext_node(model, b.NODENO));
    }
    const bool have_spc = !spc_groups.empty() || any_prescribed_flag;
    const bool have_disps = any_prescribed_flag;
    const long spc_sid = 1;

    // BLDEP -> MPC: group terms by (dependent node, dependent dof).
    // Dependent dofs flagged 3 on BNBCD but missing from BLDEP are
    // reported as a warning.
    const bool have_mpc = !model.ldeps.empty();
    const long mpc_sid = 1;
    if (have_mpc) {
        std::set<std::pair<long, long>> covered;   // (node, dof) pairs
        for (const auto &d : model.ldeps)
            for (long dd : d.DEPDOF) covered.insert({d.NODENO, dd});
        for (const auto &b : model.bcs) {
            const auto dofs = dofs_of(model, b.NODENO, b.NDOF);
            for (size_t i = 0; i < b.FIX.size() && i < dofs.size(); ++i)
                if (b.FIX[i] == fc::bnbcd::fix_key::LINDEP &&
                    !covered.count({b.NODENO, (long)dofs[i]}))
                    st.missing_ldep++;
        }
    } else {
        for (const auto &b : model.bcs)
            for (auto f : b.FIX)
                if (f == fc::bnbcd::fix_key::LINDEP) st.missing_ldep++;
    }

    // --------------------------------------------------------------
    // Executive / case control (plain text; ExtFEMIO covers the bulk)
    // --------------------------------------------------------------
    out << "$\n$ Converted from GeniE/Sesam FEM file: " << in_path << "\n"
        << "$ Converter: fem2bdf (ExtFEMIO, github.com/dnvgl/ExtFEMIO)\n$\n";
    for (const auto &t : model.user_texts) out << "$ TEXT: " << t << "\n";
    out << "SOL 101\nCEND\n"
        << "TITLE = Converted from " << in_path << "\n"
        << "ECHO = NONE\n"
        << "DISPLACEMENT(PRINT) = ALL\n"
        << "SPCFORCES(PRINT) = ALL\n"
        << "STRESS(PRINT) = ALL\n";

    if (load_cases.empty()) {
        out << "SUBCASE 1\n  LABEL = Empty\n";
        if (have_spc) out << "  SPC = " << spc_sid << "\n";
        if (have_mpc) out << "  MPC = " << mpc_sid << "\n";
    } else {
        size_t sub = 0;
        for (long llc : load_cases) {
            ++sub;
            std::string label = "FEM load case " + std::to_string(llc);
            auto it = model.load_names.find(llc);
            if (it != model.load_names.end()) label = it->second;
            out << "SUBCASE " << sub << "\n"
                << "  LABEL = " << label << "\n";
            // with prescribed displacements each case gets its own
            // SPC set carrying that case's values
            if (have_spc) out << "  SPC = " << (have_disps ? llc : spc_sid) << "\n";
            if (have_mpc) out << "  MPC = " << mpc_sid << "\n";
            out << "  LOAD = " << llc << "\n";
        }
    }
    out << "BEGIN BULK\n";

    // --------------------------------------------------------------
    // Materials: MISOSEL -> MAT1; MORSMEL -> MAT2 (anisotropic
    // membrane material: the D matrix maps 1:1 onto MAT2's G terms,
    // ALPHA1/ALPHA2 -> A1/A2; the Q direction vector is not
    // representable and dropped with a warning)
    // --------------------------------------------------------------
    out << "$\n$ Materials (MISOSEL -> MAT1, MORSMEL -> MAT2)\n$\n";
    for (const auto &kv : model.mats) {
        long mid = kv.first;
        double E = kv.second.YOUNG;
        double nu = kv.second.POISS;
        double rho = kv.second.RHO;
        // G left blank -> Nastran derives it from E and NU.
        bc::mat1 card(&mid, &E, nullptr, &nu, &rho);
        out << card;
    }
    for (const auto &kv : model.morsmels) {
        long mid = kv.first;
        const auto &m = kv.second;
        double G11 = m.D11, G12 = m.D21, G13 = m.D31;
        double G22 = m.D22, G23 = m.D32, G33 = m.D33;
        double rho = m.RHO;
        double a1 = m.ALPHA1, a2 = m.ALPHA2;
        bc::mat2 card(&mid, &G11, &G12, &G13, &G22, &G23, &G33, &rho,
                      &a1, &a2, nullptr, nullptr, nullptr);
        out << card;
        if (m.Q1 != 0.0 || m.Q2 != 0.0 || m.Q3 != 0.0) st.morsmel_q++;
    }

    // --------------------------------------------------------------
    // Nodes: GNODE/GCOORD -> GRID (external numbering)
    // --------------------------------------------------------------
    out << "$\n$ Nodes (GNODE/GCOORD -> GRID)\n$\n";
    size_t n_grid = 0;
    for (const auto &kv : model.coords) {
        long id = ext_node(model, kv.first);
        double x = kv.second.XCOORD;
        double y = kv.second.YCOORD;
        double z = kv.second.ZCOORD;
        bc::grid card(&id, nullptr, &x, &y, &z);
        out << card;
        ++n_grid;
    }

    // --------------------------------------------------------------
    // Property usage scan
    // --------------------------------------------------------------
    long max_id = 1;
    for (const auto &kv : model.thickness) max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.beamg)     max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.pipes)     max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.iorh)      max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.iorhrs)    max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.usyis)     max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.barm)      max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.boxes)     max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.dobos)     max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.chans)     max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.chanrs)    max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.lsecs)     max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.lsecrs)    max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.tonps)     max_id = std::max(max_id, kv.first);
    for (const auto &kv : model.springs)   max_id = std::max(max_id, kv.first);
    PidPool pids(max_id + 1);

    // Determine which properties are actually used, per element
    // category. Shell/beam/rod elements need a geometry reference
    // (GEONO); solids need only a material; springs/mass elements
    // reference MGSPRNG/MGMASS through the material number.
    struct Use {
        long geono, matno;
    };
    std::vector<Use> shell_use, beam_use, rod_use;
    std::set<long> solid_mats;    ///< MATNOs used by solid elements
    std::set<long> spring_props;  ///< MGSPRNG ids used by GSPR elements
    std::set<long> mass_props;    ///< MGMASS ids used by GMAS elements
    std::set<std::pair<long, long>> seen_use;
    for (const auto &kv : model.elems) {
        const auto &e = kv.second;
        auto rit = model.refs.find(e.ELNO);
        if (rit == model.refs.end()) continue;
        const long geono = elem_geono(rit->second);
        const long matno = rit->second.MATNO;
        if (matno != 0 && !model.mats.count(matno) &&
            !model.morsmels.count(matno) &&
            !model.springs.count(matno) && !model.mgmass.count(matno))
            st.missing_mat++;

        if (is_solid(e.ELTYP)) {
            solid_mats.insert(matno);
            continue;
        }
        if (e.ELTYP == fe::el_types::GSPR) {
            if (matno == 0 || !model.springs.count(matno)) {
                st.missing_spring++;
            } else {
                spring_props.insert(matno);
            }
            continue;
        }
        if (e.ELTYP == fe::el_types::GMAS) {
            if (matno == 0 || !model.mgmass.count(matno)) {
                st.missing_mass++;
            } else {
                mass_props.insert(matno);
            }
            continue;
        }
        // No geometry reference at all: PID would be 0, which Nastran
        // rejects. Count here (the element loop skips them too).
        if (geono == 0) { st.no_geono++; continue; }
        const auto key = std::make_pair(geono, matno);
        const bool fresh = seen_use.insert(key).second;
        switch (e.ELTYP) {
        case fe::el_types::FQUS_FFQ:
        case fe::el_types::FTRS_FFTR:
            if (fresh) shell_use.push_back({geono, matno});
            break;
        case fe::el_types::BEAS:
        case fe::el_types::BEPS:
            if (fresh) beam_use.push_back({geono, matno});
            break;
        case fe::el_types::TESS:
            if (fresh) rod_use.push_back({geono, matno});
            break;
        default:
            break;
        }
    }

    // Materials referenced by elements but with no MISOSEL/MORSMEL card
    // in the input (MATNO=0 means "no material assigned" in GeniE).
    // Nastran requires a valid MID, so emit an all-zero placeholder MAT1
    // per missing MATNO and warn.  Zeros are deliberate: the model
    // cannot be solved until corrected, so nothing silently plausible.
    std::map<long, long> mid_fix;   // MATNO -> MID actually written
    {
        auto has_mat = [&](long m) {
            return m != 0 &&
                (model.mats.count(m) || model.morsmels.count(m));
        };
        std::set<long> missing_mats;
        for (const auto &u : shell_use)
            if (!has_mat(u.matno)) missing_mats.insert(u.matno);
        for (const auto &u : beam_use)
            if (!has_mat(u.matno)) missing_mats.insert(u.matno);
        for (const auto &u : rod_use)
            if (!has_mat(u.matno)) missing_mats.insert(u.matno);
        for (long m : solid_mats)
            if (!has_mat(m)) missing_mats.insert(m);
        if (!missing_mats.empty()) {
            out << "$\n$ Placeholder materials (no MISOSEL in input; all-zero MAT1)\n$\n";
            for (long m : missing_mats) {
                long mid = (m != 0) ? m : max_id + 1;
                mid_fix[m] = mid;
                // All fields zero: syntactically legal, but the model
                // cannot be solved until the user fills in real values.
                // Deliberately NOT a plausible-looking default (steel),
                // so a forgotten fix can never yield believable wrong
                // results - the solver fails immediately instead.
                double E = 0.0, nu = 0.0, rho = 0.0;
                bc::mat1 card(&mid, &E, nullptr, &nu, &rho);
                out << card;
                st.default_mats++;
            }
        }
    }
    auto resolve_mid = [&](long matno) -> long {
        auto it = mid_fix.find(matno);
        return it != mid_fix.end() ? it->second : matno;
    };

    out << "$\n$ Shell properties (GELTH -> PSHELL)\n$\n";
    // GEONOs whose property card cannot be built (beams/rods without
    // any section card). Elements referencing these PIDs are skipped
    // later so the BDF never contains a dangling property reference.
    std::set<long> bad_pid;
    for (const auto &u : shell_use) {
        auto it = model.thickness.find(u.geono);
        long pid = pids.get(u.geono, u.matno);
        long mid = resolve_mid(u.matno);
        // No thickness defined for this GEONO (no GELTH card, or a
        // non-positive TH). Skipping the PSHELL would leave the shell
        // elements pointing at a non-existent PID, so emit a zero
        // thickness placeholder instead and warn loudly. t=0 is
        // syntactically legal but cannot be solved - deliberately so:
        // a forgotten fix fails immediately instead of producing
        // believable wrong results.
        double t = 0.0;
        if (it == model.thickness.end() || it->second.TH <= 0.0) {
            st.default_thk++;
        } else {
            t = it->second.TH;
        }
        bc::pshell card(&pid, &mid, &t);
        out << card;
    }

    out << "$\n$ Solid properties (material -> PSOLID)\n$\n";
    for (long m : solid_mats) {
        long pid = pids.get(m, m);
        out << "PSOLID  " << f8i(pid) << f8i(resolve_mid(m)) << "\n";
    }

    out << "$\n$ Beam properties (GBEAMG -> PBEAM; GPIPE/GIORH/GIORHR/\n"
        << "$   GUSYI/GBOX/GDOBO/GCHAN/GCHANR/GLSEC/GLSECR/GTONP/GBARM\n"
        << "$   -> PBEAML)\n$\n";
    for (const auto &u : beam_use) {
        long pid = pids.get(u.geono, u.matno);
        long mid = resolve_mid(u.matno);

        auto git = model.beamg.find(u.geono);
        auto pit = model.pipes.find(u.geono);
        auto iit = model.iorh.find(u.geono);
        auto rit = model.iorhrs.find(u.geono);
        auto uit = model.usyis.find(u.geono);
        auto bit = model.barm.find(u.geono);
        auto xit = model.boxes.find(u.geono);
        auto dit = model.dobos.find(u.geono);
        auto cit = model.chans.find(u.geono);
        auto crit = model.chanrs.find(u.geono);
        auto lit = model.lsecs.find(u.geono);
        auto lrit = model.lsecrs.find(u.geono);
        auto tit = model.tonps.find(u.geono);

        const std::string group("MSCBML0");
        // Precedence: shape-specific section cards first.  GeniE writes
        // a GBEAMG (generic computed data) alongside the shape card
        // (GLSEC, GIORH, ...) under the SAME GEONO, so GBEAMG must be
        // the fallback, not the first match.
        if (pit != model.pipes.end()) {
            // GPIPE -> PBEAML "TUBE" (DIM1 = outer diameter, DIM2 = wall)
            const auto &s = pit->second;
            double t = s.T > 0.0 ? s.T : (s.DY - s.DI) / 2.0;
            const std::string type("TUBE");
            std::vector<std::vector<double>> dim{{s.DY, t}};
            bc::pbeaml card(&pid, &mid, &group, &type, &dim);
            out << card;
        } else if (iit != model.iorh.end()) {
            // GIORH -> PBEAML "I"
            // DIM1=HZ(height) DIM2=BB(bottom flange width)
            // DIM3=BT(top flange width) DIM4=TY(web thickness)
            // DIM5=TB(bottom flange thickness) DIM6=TT(top flange thk)
            const auto &s = iit->second;
            out << "$ GIORH GEONO=" << u.geono << " -> PBEAML I\n";
            const std::string type("I");
            std::vector<std::vector<double>> dim{
                {s.HZ, s.BB, s.BT, s.TY, s.TB, s.TT}};
            bc::pbeaml card(&pid, &mid, &group, &type, &dim);
            out << card;
        } else if (rit != model.iorhrs.end()) {
            // GIORHR (I/H with inside curvature) -> PBEAML "I";
            // fillet radii RT/RB are not representable, dropped
            const auto &s = rit->second;
            st.curvature_ignored++;
            out << "$ GIORHR GEONO=" << u.geono
                << " -> PBEAML I (fillet radii dropped)\n";
            const std::string type("I");
            std::vector<std::vector<double>> dim{
                {s.HZ, s.BB, s.BT, s.TY, s.TB, s.TT}};
            bc::pbeaml card(&pid, &mid, &group, &type, &dim);
            out << card;
        } else if (uit != model.usyis.end()) {
            // GUSYI (unsymmetric I) -> PBEAML "I"; the flange
            // half-widths B1/B2 (off-centre web) are not
            // representable: keep full widths, drop the offset
            const auto &s = uit->second;
            if (std::fabs(s.B1 - 0.5 * s.BT) > 1e-9 * std::max({s.BT, 1.0}) ||
                std::fabs(s.B2 - 0.5 * s.BB) > 1e-9 * std::max({s.BB, 1.0})) {
                st.usyi_offset++;
                out << "$ WARNING: GUSYI GEONO=" << u.geono
                    << " has off-centre web; flange offsets dropped\n";
            }
            out << "$ GUSYI GEONO=" << u.geono << " -> PBEAML I\n";
            const std::string type("I");
            std::vector<std::vector<double>> dim{
                {s.HZ, s.BB, s.BT, s.TY, s.TB, s.TT}};
            bc::pbeaml card(&pid, &mid, &group, &type, &dim);
            out << card;
        } else if (dit != model.dobos.end()) {
            // GDOBO (double bottom: two plates + web) -> PBEAML "I"
            // with equal flanges of width BY
            const auto &s = dit->second;
            out << "$ GDOBO GEONO=" << u.geono << " -> PBEAML I\n";
            const std::string type("I");
            std::vector<std::vector<double>> dim{
                {s.HZ, s.BY, s.BY, s.TY, s.TB, s.TT}};
            bc::pbeaml card(&pid, &mid, &group, &type, &dim);
            out << card;
        } else if (tit != model.tonps.end()) {
            // GTONP (T on plate) -> PBEAML "I": the effective plate
            // (BP x TP) acts as the bottom flange, the T flange
            // (BT x TT) as the top flange.  HZ is taken as the
            // height of the T part (plate excluded).
            const auto &s = tit->second;
            st.tonp_as_i++;
            out << "$ GTONP GEONO=" << u.geono
                << " -> PBEAML I (plate as bottom flange)\n";
            const std::string type("I");
            std::vector<std::vector<double>> dim{
                {s.HZ + s.TP, s.BP, s.BT, s.TY, s.TP, s.TT}};
            bc::pbeaml card(&pid, &mid, &group, &type, &dim);
            out << card;
        } else if (cit != model.chans.end() || crit != model.chanrs.end()) {
            // GCHAN / GCHANR -> PBEAML "CHAN"
            // DIM1=flange width DIM2=height DIM3=flange thk DIM4=web thk
            double hz, ty, by, tz;
            const char *cname;
            if (cit != model.chans.end()) {
                hz = cit->second.HZ; ty = cit->second.TY;
                by = cit->second.BY; tz = cit->second.TZ;
                cname = "GCHAN";
            } else {
                hz = crit->second.HZ; ty = crit->second.TY;
                by = crit->second.BY; tz = crit->second.TZ;
                cname = "GCHANR";
                st.curvature_ignored++;
            }
            out << "$ " << cname << " GEONO=" << u.geono << " -> PBEAML CHAN"
                << (crit != model.chanrs.end() && cit == model.chans.end()
                    ? " (fillet radius dropped)" : "") << "\n";
            const std::string type("CHAN");
            std::vector<std::vector<double>> dim{{by, hz, tz, ty}};
            bc::pbeaml card(&pid, &mid, &group, &type, &dim);
            out << card;
        } else if (lit != model.lsecs.end() || lrit != model.lsecrs.end()) {
            // GLSEC / GLSECR -> PBEAML "L"
            // DIM1=horizontal leg width DIM2=vertical leg height
            // DIM3=horizontal leg thk DIM4=vertical leg thk
            double hz, ty, by, tz;
            const char *cname;
            if (lit != model.lsecs.end()) {
                hz = lit->second.HZ; ty = lit->second.TY;
                by = lit->second.BY; tz = lit->second.TZ;
                cname = "GLSEC";
            } else {
                hz = lrit->second.HZ; ty = lrit->second.TY;
                by = lrit->second.BY; tz = lrit->second.TZ;
                cname = "GLSECR";
                st.curvature_ignored++;
            }
            out << "$ " << cname << " GEONO=" << u.geono << " -> PBEAML L"
                << (lrit != model.lsecrs.end() && lit == model.lsecs.end()
                    ? " (fillet radius dropped)" : "") << "\n";
            const std::string type("L");
            std::vector<std::vector<double>> dim{{by, hz, tz, ty}};
            bc::pbeaml card(&pid, &mid, &group, &type, &dim);
            out << card;
        } else if (xit != model.boxes.end()) {
            // GBOX -> PBEAML "BOX" (equal flanges) or "BOX1"
            const auto &s = xit->second;
            std::vector<std::vector<double>> dim;
            std::string type;
            if (std::fabs(s.TT - s.TB) <= 1e-9 * std::max({s.TT, s.TB, 1.0})) {
                // BOX: DIM1=width DIM2=height DIM3=flange thk DIM4=web thk
                type = "BOX";
                dim.push_back({s.BY, s.HZ, s.TT, s.TY});
            } else {
                // BOX1: DIM3=top thk DIM4=bottom thk DIM5/DIM6=web thk
                type = "BOX1";
                dim.push_back({s.BY, s.HZ, s.TT, s.TB, s.TY, s.TY});
            }
            out << "$ GBOX GEONO=" << u.geono << " -> PBEAML " << type << "\n";
            bc::pbeaml card(&pid, &mid, &group, &type, &dim);
            out << card;
        } else if (bit != model.barm.end()) {
            const auto &s = bit->second;
            if (std::fabs(s.BT - s.BB) <= 1e-9 * std::max({s.BT, s.BB, 1.0})) {
                // rectangular bar -> PBEAML "BAR" (DIM1=width, DIM2=height)
                out << "$ GBARM GEONO=" << u.geono << " -> PBEAML BAR\n";
                const std::string type("BAR");
                std::vector<std::vector<double>> dim{{s.BB, s.HZ}};
                bc::pbeaml card(&pid, &mid, &group, &type, &dim);
                out << card;
            } else {
                // trapezoidal bar has no native Nastran shape:
                // fall back to equivalent PBEAM
                st.trapezoid_bar++;
                out << "$ GBARM GEONO=" << u.geono
                    << " trapezoidal; converted to equivalent PBEAM\n";
                SectionProps sp = props_gbarm(s);
                std::vector<double> A{sp.area}, I1{sp.iy}, I2{sp.iz},
                    I12{0.0}, J{sp.j};
                bc::pbeam card(&pid, &mid, &A, &I1, &I2, &I12, &J);
                out << card;
            }
        } else if (git != model.beamg.end()) {
            // GBEAMG -> PBEAM (direct property copy); fallback used
            // only when no shape-specific card exists for this GEONO.
            const auto &s = git->second;
            std::vector<double> A{s.AREA}, I1{s.IY}, I2{s.IZ}, J{s.IX};
            // blank I12 when no product of inertia (keeps small-field
            // output and matches the Nastran default of 0)
            std::vector<double> I12;
            std::vector<double> *pI12 = nullptr;
            if (s.IYZ != 0.0) {
                I12.push_back(-s.IYZ);
                pI12 = &I12;
            }
            double k1 = 0.0, k2 = 0.0;
            double *pk1 = nullptr, *pk2 = nullptr;
            if (s.AREA > 0.0 && s.SHARZ > 0.0) { k1 = s.SHARZ / s.AREA; pk1 = &k1; }
            if (s.AREA > 0.0 && s.SHARY > 0.0) { k2 = s.SHARY / s.AREA; pk2 = &k2; }
            bc::pbeam card(&pid, &mid, &A, &I1, &I2, pI12, &J,
                           nullptr, nullptr, nullptr, nullptr, nullptr,
                           nullptr, nullptr, nullptr, nullptr, nullptr,
                           nullptr, pk1, pk2);
            out << card;
        } else {
            st.missing_geono++;
            bad_pid.insert(pid);
            warn_skipped(st, "beam GEONO " + std::to_string(u.geono));
            continue;
        }
    }

    out << "$\n$ Rod properties (GBEAMG/GPIPE -> PROD)\n$\n";
    for (const auto &u : rod_use) {
        long pid = pids.get(u.geono, u.matno);
        long mid = resolve_mid(u.matno);
        double area = 0.0;
        auto git = model.beamg.find(u.geono);
        auto pit = model.pipes.find(u.geono);
        if (git != model.beamg.end()) area = git->second.AREA;
        else if (pit != model.pipes.end()) area = pipe_area(pit->second);
        else {
            st.missing_geono++;
            bad_pid.insert(pid);
            warn_skipped(st, "rod GEONO " + std::to_string(u.geono));
            continue;
        }
        bc::prod card(&pid, &mid, &area);
        out << card;
    }

    // Spring properties: MGSPRNG -> PBUSH (diagonal terms only)
    out << "$\n$ Spring properties (MGSPRNG -> PBUSH)\n$\n";
    for (long key : spring_props) {
        const auto &s = model.springs[key];
        long pid = pids.get(key, 0);
        std::vector<double> k;
        double max_diag = 0.0, max_off = 0.0;
        for (size_t i = 0; i < s.K.size(); ++i) {
            for (size_t j = 0; j < s.K[i].size(); ++j) {
                if (i == j)
                    max_diag = std::max(max_diag, std::fabs(s.K[i][j]));
                else
                    max_off = std::max(max_off, std::fabs(s.K[i][j]));
            }
        }
        for (size_t i = 0; i < 6; ++i)
            k.push_back(i < s.K.size() && i < s.K[i].size() ? s.K[i][i] : 0.0);
        if (max_off > 1e-6 * std::max(1.0, max_diag)) {
            st.spring_coupling++;
            out << "$ WARNING: MGSPRNG MATNO=" << key
                << " has off-diagonal stiffness; coupling terms dropped\n";
        }
        bc::pbush card(&pid, &k);
        out << card;
    }

    // --------------------------------------------------------------
    // Elements
    // --------------------------------------------------------------
    out << "$\n$ Elements\n$\n";
    size_t n_quad = 0, n_tria = 0, n_beam = 0, n_rod = 0;
    size_t n_hex = 0, n_tet = 0, n_pri = 0, n_bush = 0, n_conm1 = 0;
    long max_eid = 0;
    bool warned_default_orient = false;

    for (const auto &kv : model.elems) {
        const auto &e = kv.second;
        auto rit = model.refs.find(e.ELNO);
        if (rit == model.refs.end()) {
            st.skipped_eltypes[static_cast<long>(e.ELTYP)]++;
            continue;
        }
        const auto &r = rit->second;
        const long geono = elem_geono(r);
        const long matno = r.MATNO;

        long eid = e.ELNOX > 0 ? e.ELNOX : e.ELNO;
        max_eid = std::max(max_eid, eid);
        long pid = 0;
        if (!is_solid(e.ELTYP) && e.ELTYP != fe::el_types::GSPR &&
            e.ELTYP != fe::el_types::GMAS) {
            if (geono == 0) continue;  // already counted in the scan pass
            pid = pids.get(geono, matno);
            // Property could not be built (e.g. beam GEONO without any
            // section card): skip the element rather than writing a
            // dangling PID reference that Nastran would reject.
            if (bad_pid.count(pid)) { st.no_prop++; continue; }
        }

        switch (e.ELTYP) {
        case fe::el_types::FQUS_FFQ: {
            if (e.NODIN.size() < 4) { st.skipped_eltypes[24]++; break; }
            long g1 = ext_node(model, e.NODIN[0]);
            long g2 = ext_node(model, e.NODIN[1]);
            long g3 = ext_node(model, e.NODIN[2]);
            long g4 = ext_node(model, e.NODIN[3]);
            bc::cquad4 card(&eid, &pid, &g1, &g2, &g3, &g4);
            out << card;
            ++n_quad;
            break;
        }
        case fe::el_types::FTRS_FFTR: {
            if (e.NODIN.size() < 3) { st.skipped_eltypes[25]++; break; }
            long g1 = ext_node(model, e.NODIN[0]);
            long g2 = ext_node(model, e.NODIN[1]);
            long g3 = ext_node(model, e.NODIN[2]);
            bc::ctria3 card(&eid, &pid, &g1, &g2, &g3);
            out << card;
            ++n_tria;
            break;
        }
        case fe::el_types::BEAS:
        case fe::el_types::BEPS: {
            if (e.NODIN.size() < 2) { st.skipped_eltypes[15]++; break; }
            long ga = ext_node(model, e.NODIN[0]);
            long gb = ext_node(model, e.NODIN[1]);

            // Orientation: the Nastran orientation vector describes the
            // element local y axis, which coincides with the Sesam local
            // z axis given by GUNIVEC.
            double x1 = 0.0, x2 = 0.0, x3 = 1.0;
            const long transno = elem_transno(r);
            auto tit = model.trans.find(transno);
            if (tit != model.trans.end()) {
                x1 = tit->second.UNIX;
                x2 = tit->second.UNIY;
                x3 = tit->second.UNIZ;
            } else if (!warned_default_orient) {
                warned_default_orient = true;
                std::cerr << "WARNING: beam(s) without GUNIVEC "
                    "orientation; using default (0,0,1) with fallback\n";
            }
            if (tit == model.trans.end()) {
                // avoid orientation vector parallel to the beam axis
                const auto &ca = model.coords[e.NODIN[0]];
                const auto &cb = model.coords[e.NODIN[1]];
                const double dx = cb.XCOORD - ca.XCOORD;
                const double dy = cb.YCOORD - ca.YCOORD;
                const double dz = cb.ZCOORD - ca.ZCOORD;
                const double len = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (len > 0.0 && std::abs(dz / len) > 0.999) {
                    x1 = 0.0; x2 = 1.0; x3 = 0.0;
                }
            }

            // Element local axes in global coordinates: x along the
            // beam (node A -> node B), z = GUNIVEC component normal to
            // x, y = z cross x.  Needed for GECC (local eccentricity)
            // rotation and for the shear-centre shift below.
            double ux = 0.0, uy = 0.0, uz = 0.0;
            double yx = 0.0, yy = 0.0, yz = 0.0;
            double zx = 0.0, zy = 0.0, zz = 1.0;
            bool axes_ok = false;
            {
                const auto &ca = model.coords[e.NODIN[0]];
                const auto &cb = model.coords[e.NODIN[1]];
                double dx = cb.XCOORD - ca.XCOORD;
                double dy = cb.YCOORD - ca.YCOORD;
                double dz = cb.ZCOORD - ca.ZCOORD;
                const double len = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (len > 0.0) {
                    ux = dx / len; uy = dy / len; uz = dz / len;
                    double tx = x1, ty = x2, tz = x3;
                    const double dot = tx * ux + ty * uy + tz * uz;
                    tx -= dot * ux; ty -= dot * uy; tz -= dot * uz;
                    const double tl = std::sqrt(tx * tx + ty * ty + tz * tz);
                    if (tl > 1e-12) {
                        zx = tx / tl; zy = ty / tl; zz = tz / tl;
                        yx = zy * uz - zz * uy;
                        yy = zz * ux - zx * uz;
                        yz = zx * uy - zy * ux;
                        axes_ok = true;
                    }
                }
            }

            // GECCEN / GECC -> CBEAM end offsets WA / WB.  GECCEN is
            // given in the global (superelement) system; GECC is given
            // in the element local system and rotated to global here.
            // Both point from the mesh node towards the element node;
            // Sesam places the element nodes on the neutral-axis line.
            auto ecc_at = [&](size_t iend, double &ex, double &ey,
                              double &ez) -> bool {
                ex = ey = ez = 0.0;
                long eccno = 0;
                if (r.ECCNO_OPT > 0) eccno = r.ECCNO_OPT;
                else if (r.ECCNO_OPT == -1 && r.ECCNO.size() > iend)
                    eccno = r.ECCNO[iend];
                if (eccno == 0) return false;
                auto gi = model.eccs.find(eccno);
                if (gi != model.eccs.end()) {
                    ex = gi->second.EX; ey = gi->second.EY; ez = gi->second.EZ;
                } else {
                    auto li = model.leccs.find(eccno);
                    if (li == model.leccs.end() || !axes_ok) return false;
                    ex = li->second.EX * ux + li->second.EY * yx +
                         li->second.EZ * zx;
                    ey = li->second.EX * uy + li->second.EY * yy +
                         li->second.EZ * zy;
                    ez = li->second.EX * uz + li->second.EY * yz +
                         li->second.EZ * zz;
                    st.gecc_local++;
                }
                return (ex != 0.0 || ey != 0.0 || ez != 0.0);
            };
            double exa, eya, eza, exb, eyb, ezb;
            const bool has_a = ecc_at(0, exa, eya, eza);
            const bool has_b = ecc_at(1, exb, eyb, ezb);
            // single eccentricity for the whole element -> both ends
            double w1a = exa, w2a = eya, w3a = eza;
            double w1b = has_b ? exb : exa;
            double w2b = has_b ? eyb : eya;
            double w3b = has_b ? ezb : eza;
            bool has_off = has_a || has_b || r.ECCNO_OPT > 0;

            // BELFIX -> CBEAM pin flags PA/PB.  FIXATION: A(i)=0 fully
            // released, =1 fully connected; fractional values and
            // OPT=SPRING rates cannot be represented by pin flags.
            std::vector<int> pa, pb;
            auto hinge_at = [&](size_t iend, std::vector<int> &pins) {
                long fixno = 0;
                if (r.FIXNO_OPT > 0) fixno = r.FIXNO_OPT;
                else if (r.FIXNO_OPT == -1 && r.FIXNO.size() > iend)
                    fixno = r.FIXNO[iend];
                if (fixno == 0) return;
                auto it = model.hinges.find(fixno);
                if (it == model.hinges.end()) { st.missing_hinge++; return; }
                const auto &h = it->second;
                if (h.TRANO != 0) st.hinge_local++;
                const bool spring =
                    (h.OPT == fc::belfix::n_opt::SPRING ||
                     h.OPT == fc::belfix::n_opt::SPRING_END);
                bool counted_spring = false;
                for (size_t i = 0; i < h.A.size() && i < 6; ++i) {
                    const double a = h.A[i];
                    if (spring) {
                        // Ci = -1: rigid; Ci = 0: free; Ci > 0: elastic
                        if (a == -1.0) continue;
                        if (!counted_spring) { st.hinge_spring++; counted_spring = true; }
                        if (a == 0.0) pins.push_back((int)i + 1);
                        else st.hinge_partial++;
                    } else {
                        if (a >= 1.0) continue;         // fully connected
                        if (a == 0.0) pins.push_back((int)i + 1);
                        else st.hinge_partial++;        // 0 < a < 1
                    }
                }
                if (!pins.empty()) st.hinge_applied++;
            };
            hinge_at(0, pa);
            hinge_at(1, pb);

            // Nastran draws PBEAML shapes from the SHEAR CENTRE (the
            // CBEAM element axis is the shear-centre axis), while Sesam
            // puts the element nodes on the centroid.  Shift the offsets
            // by the centroid->shear-centre vector (GBEAMG SHCENY/SHCENZ,
            // given in the Sesam local y/z axes) rotated to global.
            {
                const long geono = elem_geono(r);
                auto git = model.beamg.find(geono);
                if (git != model.beamg.end() && shape_card_for(model, geono) &&
                    axes_ok &&
                    (git->second.SHCENY != 0.0 || git->second.SHCENZ != 0.0)) {
                    const double oy = git->second.SHCENY;
                    const double oz = git->second.SHCENZ;
                    const double gx = oy * yx + oz * zx;
                    const double gy = oy * yy + oz * zy;
                    const double gz = oy * yz + oz * zz;
                    w1a += gx; w2a += gy; w3a += gz;
                    w1b += gx; w2b += gy; w3b += gz;
                    has_off = true;
                    st.shcen_adjust++;
                }
            }

            if (has_off || !pa.empty() || !pb.empty()) {
                bc::cbeam card(&eid, &pid, &ga, &gb, &x1, &x2, &x3,
                               (std::string*)nullptr,
                               pa.empty() ? nullptr : &pa,
                               pb.empty() ? nullptr : &pb,
                               has_off ? &w1a : nullptr,
                               has_off ? &w2a : nullptr,
                               has_off ? &w3a : nullptr,
                               has_off ? &w1b : nullptr,
                               has_off ? &w2b : nullptr,
                               has_off ? &w3b : nullptr);
                out << card;
            } else {
                bc::cbeam card(&eid, &pid, &ga, &gb, &x1, &x2, &x3);
                out << card;
            }
            ++n_beam;
            break;
        }
        case fe::el_types::TESS: {
            if (e.NODIN.size() < 2) { st.skipped_eltypes[10]++; break; }
            long g1 = ext_node(model, e.NODIN[0]);
            long g2 = ext_node(model, e.NODIN[1]);
            bc::crod card(&eid, &pid, &g1, &g2);
            out << card;
            ++n_rod;
            break;
        }
        case fe::el_types::LHEX:
        case fe::el_types::IHEX:
        case fe::el_types::TETR:
        case fe::el_types::ITET:
        case fe::el_types::TPRI:
        case fe::el_types::IPRI: {
            const bool is_hex = (e.ELTYP == fe::el_types::LHEX ||
                                 e.ELTYP == fe::el_types::IHEX);
            const bool is_tet = (e.ELTYP == fe::el_types::TETR ||
                                 e.ELTYP == fe::el_types::ITET);
            const size_t nc = is_hex ? 8 : (is_tet ? 4 : 6);
            if (e.NODIN.size() < nc) {
                st.skipped_eltypes[static_cast<long>(e.ELTYP)]++;
                break;
            }
            pid = pids.get(matno, matno);   // PSOLID id

            // external ids of the corner nodes in Sesam local order
            std::vector<long> cn;
            for (size_t i = 0; i < nc; ++i)
                cn.push_back(ext_node(model, e.NODIN[i]));

            // orientation check on the corner geometry
            const auto cc = corner_coords(model, e.NODIN, nc);
            if (!cc.empty()) {
                const double d = is_hex ? hex_detJ(cc)
                               : is_tet ? tet_vol6(cc)
                                        : wedge_detJ(cc);
                if (d < 0.0) {
                    if (is_hex) { std::swap(cn[1], cn[3]);
                                  std::swap(cn[5], cn[7]); }
                    else if (is_tet) { std::swap(cn[1], cn[2]); }
                    else { std::swap(cn[1], cn[2]);
                           std::swap(cn[4], cn[5]); }
                }
            }

            // mid-edge nodes (quadratic elements): Sesam and Nastran
            // number them differently; resolve by edge endpoints.
            std::vector<long> outnodes = cn;
            const size_t nn = e.NODIN.size();
            if (nn > nc) {
                // unordered corner pair -> mid node (Sesam layout)
                auto edge_key = [](long a, long b) {
                    if (a > b) std::swap(a, b);
                    return std::make_pair(a, b);
                };
                std::map<std::pair<long, long>, long> mid_of;
                auto add_mids = [&](std::initializer_list<std::pair<int, int>> edges,
                                    size_t first) {
                    size_t k = 0;
                    for (const auto &ed : edges) {
                        const long a = ext_node(model, e.NODIN[ed.first - 1]);
                        const long b = ext_node(model, e.NODIN[ed.second - 1]);
                        if (first + k < nn)
                            mid_of[edge_key(a, b)] =
                                ext_node(model, e.NODIN[first + k]);
                        ++k;
                    }
                };
                // Nastran mid-edge order, as corner-position pairs
                std::vector<std::pair<int, int>> nas_edges;
                if (is_hex) {
                    // Sesam IHEX: 9-12 bottom, 13-16 vertical, 17-20 top
                    add_mids({{1,2},{2,3},{3,4},{4,1}}, 8);
                    add_mids({{1,5},{2,6},{3,7},{4,8}}, 12);
                    add_mids({{5,6},{6,7},{7,8},{8,5}}, 16);
                    nas_edges = {{0,1},{1,2},{2,3},{3,0},
                                 {4,5},{5,6},{6,7},{7,4},
                                 {0,4},{1,5},{2,6},{3,7}};
                } else if (is_tet) {
                    // Sesam ITET == Nastran CTETRA10 edge order
                    add_mids({{1,2},{2,3},{3,1},{1,4},{2,4},{3,4}}, 4);
                    nas_edges = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
                } else {
                    // Sesam IPRI: 7-9 bottom, 10-12 vertical, 13-15 top
                    add_mids({{1,2},{2,3},{3,1}}, 6);
                    add_mids({{1,4},{2,5},{3,6}}, 9);
                    add_mids({{4,5},{5,6},{6,4}}, 12);
                    nas_edges = {{0,1},{1,2},{2,0},
                                 {3,4},{4,5},{5,3},
                                 {0,3},{1,4},{2,5}};
                }
                bool ok = true;
                for (const auto &ed : nas_edges) {
                    auto mit = mid_of.find(edge_key(cn[ed.first], cn[ed.second]));
                    if (mit == mid_of.end()) { ok = false; break; }
                    outnodes.push_back(mit->second);
                }
                if (!ok) {
                    st.skipped_eltypes[static_cast<long>(e.ELTYP)]++;
                    break;
                }
            }

            if (is_hex) {
                out << "CHEXA   " << f8i(eid) << f8i(pid);
                for (size_t i = 0; i < outnodes.size(); ++i) {
                    if (i == 6 || i == 14) out << "\n        ";
                    out << f8i(outnodes[i]);
                }
                out << "\n";
                ++n_hex;
            } else if (is_tet) {
                out << "CTETRA  " << f8i(eid) << f8i(pid);
                for (size_t i = 0; i < outnodes.size(); ++i) {
                    if (i == 6) out << "\n        ";
                    out << f8i(outnodes[i]);
                }
                out << "\n";
                ++n_tet;
            } else {
                out << "CPENTA  " << f8i(eid) << f8i(pid);
                for (size_t i = 0; i < outnodes.size(); ++i) {
                    if (i == 6 || i == 14) out << "\n        ";
                    out << f8i(outnodes[i]);
                }
                out << "\n";
                ++n_pri;
            }
            break;
        }
        case fe::el_types::GSPR: {
            // spring to ground: MGSPRNG referenced via MATNO
            if (matno == 0 || !model.springs.count(matno)) break; // counted
            if (e.NODIN.empty()) { st.skipped_eltypes[18]++; break; }
            pid = pids.get(matno, 0);
            long ga = ext_node(model, e.NODIN[0]);
            bc::cbush card(&eid, &pid, &ga);
            out << card;
            ++n_bush;
            if (elem_transno(r) != 0) st.spring_local++;
            break;
        }
        case fe::el_types::GMAS: {
            // 1-noded mass element: MGMASS referenced via MATNO
            if (matno == 0 || !model.mgmass.count(matno)) break; // counted
            if (e.NODIN.empty()) { st.skipped_eltypes[11]++; break; }
            const auto &mm = model.mgmass[matno];
            const long ndof =
                (long)((std::sqrt(8.0 * mm.size() + 1.0) - 1.0) / 2.0);
            const auto dofs = dofs_of(model, e.NODIN[0], ndof);
            double M[6][6] = {};
            // MGMASS stores the lower triangle column by column:
            // M(1,1), M(2,1), ..., M(n,1), M(2,2), ..., M(n,n)
            size_t idx = 0;
            for (long j = 0; j < ndof && idx < mm.size(); ++j) {
                for (long i = j; i < ndof && idx < mm.size(); ++i) {
                    const int di = i < (long)dofs.size() ? dofs[i] - 1 : -1;
                    const int dj = j < (long)dofs.size() ? dofs[j] - 1 : -1;
                    if (di >= 0 && di < 6 && dj >= 0 && dj < 6) {
                        M[di][dj] = mm[idx];
                        M[dj][di] = mm[idx];
                    }
                    ++idx;
                }
            }
            // CONM1 wants the 21 lower-triangle values row by row:
            // M11, M21, M22, M31, ..., M66 (individual arguments)
            double mij[21];
            {
                int k = 0;
                for (int i = 0; i < 6; ++i)
                    for (int j = 0; j <= i; ++j)
                        mij[k++] = M[i][j];
            }
            long g = ext_node(model, e.NODIN[0]);
            long cid = 0;   // basic coordinate system
            bc::conm1 card(&eid, &g, &cid,
                           &mij[0],
                           &mij[1], &mij[2],
                           &mij[3], &mij[4], &mij[5],
                           &mij[6], &mij[7], &mij[8], &mij[9],
                           &mij[10], &mij[11], &mij[12], &mij[13], &mij[14],
                           &mij[15], &mij[16], &mij[17], &mij[18], &mij[19],
                           &mij[20]);
            out << card;
            ++n_conm1;
            if (elem_transno(r) != 0) st.spring_local++;
            break;
        }
        default:
            st.skipped_eltypes[static_cast<long>(e.ELTYP)]++;
            break;
        }
    }

    // Nodal point masses: BNMASS -> CMASS2 (one card per dof)
    size_t n_cmass2 = 0;
    if (!model.nodal_masses.empty()) {
        out << "$\n$ Nodal point masses (BNMASS -> CMASS2)\n$\n";
        for (const auto &nm : model.nodal_masses) {
            long g = ext_node(model, nm.nodeno);
            const auto dofs =
                dofs_of(model, nm.nodeno, (long)nm.mass.size());
            for (size_t i = 0; i < nm.mass.size() && i < dofs.size(); ++i) {
                if (nm.mass[i] == 0.0) continue;
                long eid = ++max_eid;
                double mv = nm.mass[i];
                std::vector<int> c1{dofs[i]};
                bc::cmass2 card(&eid, &mv, &g, &c1);
                out << card;
                ++n_cmass2;
            }
        }
    }

    // --------------------------------------------------------------
    // Boundary conditions: BNBCD -> SPC1 (fixed at zero) and
    // BNBCD/BNDISPL -> SPC (prescribed values, per load case).
    // Written directly; ExtFEMIO implements neither card.
    // --------------------------------------------------------------
    size_t n_spc = 0;
    if (have_spc) {
        out << "$\n$ Boundary conditions (BNBCD"
            << (have_disps ? " + BNDISPL -> SPC1 / SPC" : " -> SPC1")
            << ")\n$\n";

        // helper: fixed-format SPC1 for one set id
        auto write_spc1 = [&out](long sid, const std::string &comp,
                                 const std::vector<long> &nodes) {
            char head[40];
            std::snprintf(head, sizeof(head), "%-8s%8d%-8s",
                          "SPC1", (int)sid, comp.c_str());
            out << head;
            size_t col = 0;
            for (size_t i = 0; i < nodes.size(); ++i) {
                if ((i < 6 && col == 6) || (i >= 6 && (i - 6) % 8 == 0))
                    out << (i == 0 ? "" : "\n        ");
                char fld[16];
                std::snprintf(fld, sizeof(fld), "%8ld", nodes[i]);
                out << fld;
                ++col;
            }
            out << "\n";
        };

        std::vector<long> spc_sets;
        if (have_disps && !load_cases.empty())
            spc_sets.assign(load_cases.begin(), load_cases.end());
        else
            spc_sets.push_back(spc_sid);

        // prescribed values per (load case, internal node)
        std::map<std::pair<long, long>, const fc::bndispl *> dmap;
        for (const auto &d : model.disps)
            dmap[std::make_pair(d.LLC, d.NODENO)] = &d;
        std::set<const fc::bndispl *> complx_noted;

        for (long sid : spc_sets) {
            for (const auto &kv : spc_groups)
                write_spc1(sid, kv.first, kv.second);

            if (!have_disps) continue;
            for (const auto &b : model.bcs) {
                const auto dofs = dofs_of(model, b.NODENO, b.NDOF);
                for (size_t i = 0; i < b.FIX.size() && i < dofs.size(); ++i) {
                    if (b.FIX[i] != fc::bnbcd::fix_key::PRESCRIBED) continue;
                    double val = 0.0;
                    auto it = dmap.find(std::make_pair(sid, b.NODENO));
                    if (it == dmap.end()) {
                        st.missing_prescribed++;
                    } else {
                        const auto *d = it->second;
                        if (d->DTYPE != 1) {
                            st.accel_disp++;
                        } else {
                            if (i < d->RDISP.size()) val = d->RDISP[i];
                            if (d->COMPLX &&
                                complx_noted.insert(d).second)
                                st.complx_disp++;
                        }
                    }
                    out << "SPC     " << f8i(sid)
                        << f8i(ext_node(model, b.NODENO))
                        << f8i(dofs[i]) << f8r(val) << "\n";
                    ++n_spc;
                }
            }
        }
    }

    // --------------------------------------------------------------
    // Linear dependences: BLDEP -> MPC (one MPC per dependent dof;
    // u_dep - sum b_i * u_indep_i = 0, first term is the dependent
    // one per Nastran convention).  Written directly; ExtFEMIO has
    // no MPC card.
    // --------------------------------------------------------------
    if (have_mpc) {
        out << "$\n$ Linear dependences (BLDEP -> MPC)\n$\n";
        std::map<std::pair<long, long>,
                 std::vector<std::tuple<long, long, double>>> terms;
        for (const auto &d : model.ldeps) {
            const long gdep = ext_node(model, d.NODENO);
            const long gind = ext_node(model, d.CNOD);
            const size_t n = std::min({d.DEPDOF.size(), d.INDEPDOF.size(),
                                       d.b.size()});
            for (size_t i = 0; i < n; ++i)
                terms[{gdep, d.DEPDOF[i]}].emplace_back(
                    gind, d.INDEPDOF[i], d.b[i]);
        }
        for (const auto &kv : terms) {
            out << "MPC     " << f8i(mpc_sid)
                << f8i(kv.first.first) << f8i(kv.first.second) << f8r(1.0);
            size_t col = 1;   // terms already on this line
            for (const auto &t : kv.second) {
                if (col == 2) { out << "\n        "; col = 0; }
                out << f8i(std::get<0>(t)) << f8i(std::get<1>(t))
                    << f8r(-std::get<2>(t));
                ++col;
            }
            out << "\n";
            ++st.n_mpc;
        }
    }

    // --------------------------------------------------------------
    // Nodal loads: BNLOAD -> FORCE / MOMENT
    // --------------------------------------------------------------
    if (!model.loads.empty()) {
        out << "$\n$ Nodal loads (BNLOAD -> FORCE / MOMENT)\n$\n";
        for (const auto &l : model.loads) {
            long sid = l.LLC;
            long g = ext_node(model, l.NODENO);
            const auto dofs = dofs_of(model, l.NODENO, l.NDOF);
            double f[3] = {0, 0, 0}, mo[3] = {0, 0, 0};
            for (size_t i = 0; i < l.RLOAD.size() && i < dofs.size(); ++i) {
                const int d = dofs[i];
                if (d >= 1 && d <= 3) f[d - 1] = l.RLOAD[i];
                else if (d >= 4 && d <= 6) mo[d - 4] = l.RLOAD[i];
            }
            const double one = 1.0;
            long cid = 0; // basic coordinate system
            if (f[0] != 0.0 || f[1] != 0.0 || f[2] != 0.0) {
                bc::force card(&sid, &g, &cid, &one, &f[0], &f[1], &f[2]);
                out << card;
            }
            if (mo[0] != 0.0 || mo[1] != 0.0 || mo[2] != 0.0) {
                bc::moment card(&sid, &g, &cid, &one, &mo[0], &mo[1], &mo[2]);
                out << card;
            }
        }
    }

    // --------------------------------------------------------------
    // Surface loads: BEUSLO (normal pressure) -> PLOAD4
    // --------------------------------------------------------------
    size_t n_pload4 = 0;
    if (!model.sloads.empty()) {
        out << "$\n$ Surface loads (BEUSLO -> PLOAD4)\n$\n";

        auto write_pload4 = [&out](long sid, long eid,
                                   const std::vector<double> &p,
                                   long g1, long g4) {
            bool uniform = true;
            for (size_t i = 1; i < p.size(); ++i)
                if (p[i] != p[0]) { uniform = false; break; }
            out << "PLOAD4  " << f8i(sid) << f8i(eid);
            if (uniform) {
                out << f8r(p.empty() ? 0.0 : p[0]);
                if (g1 || g4) out << "                        "; // P2..P4
            } else {
                for (size_t i = 0; i < 4; ++i)
                    out << (i < p.size() ? f8r(p[i])
                                         : std::string("        "));
            }
            if (g1 || g4)
                out << f8i(g1) << (g4 ? f8i(g4) : std::string("        "));
            out << "\n";
        };

        for (const auto &l : model.sloads) {
            const long lotyp = std::labs(l.LOTYP);
            if (lotyp == 3) { st.hydro_dummy++; continue; }
            if (lotyp == 2) { st.comp_form_load++; continue; }
            if (lotyp != 1) { st.bad_side++; continue; }
            if (l.LOTYP < 0) st.follower_load++;
            if (l.COMPLX) st.complx_load++;

            auto eit = model.elems.find(l.ELNO);
            if (eit == model.elems.end()) { st.missing_elem_load++; continue; }
            const auto &e = eit->second;
            const long eid = e.ELNOX > 0 ? e.ELNOX : e.ELNO;
            const long sid = l.LLC;

            if (e.ELTYP == fe::el_types::FQUS_FFQ ||
                e.ELTYP == fe::el_types::FTRS_FFTR) {
                // shell: SIDE 1 = bottom face (+p), 3 = top face (-p),
                // 2 = mid-surface (treated like SIDE 3)
                double sign = 1.0;
                if (l.SIDE == 2) { st.side_mid++; sign = -1.0; }
                else if (l.SIDE != 1) sign = -1.0;
                const size_t nn =
                    e.ELTYP == fe::el_types::FQUS_FFQ ? 4 : 3;
                std::vector<double> p;
                for (size_t i = 0; i < nn; ++i)
                    p.push_back(sign * (i < l.RLOADi.size()
                                        ? l.RLOADi[i] : 0.0));
                write_pload4(sid, eid, p, 0, 0);
                ++n_pload4;
            } else if (is_solid(e.ELTYP)) {
                const bool is_hex = (e.ELTYP == fe::el_types::LHEX ||
                                     e.ELTYP == fe::el_types::IHEX);
                const bool is_tet = (e.ELTYP == fe::el_types::TETR ||
                                     e.ELTYP == fe::el_types::ITET);
                const size_t nc = is_hex ? 8 : (is_tet ? 4 : 6);
                if (e.NODIN.size() < nc) { st.missing_elem_load++; continue; }
                const int ndig = is_hex ? 6 : (is_tet ? 4 : 5);
                // SIDE is a digit flag: the digit position from the
                // left (padded to ndig) gives the loaded side number
                if (l.SIDE <= 0) { st.bad_side++; continue; }
                long p10 = 1;
                int ndigits = 0;
                while (p10 <= l.SIDE) { p10 *= 10; ++ndigits; }
                if (p10 / 10 != l.SIDE) { st.bad_side++; continue; }
                const int side = ndig - ndigits + 1;
                const int nface = is_hex ? 6 : (is_tet ? 4 : 5);
                if (side < 1 || side > nface) { st.bad_side++; continue; }

                std::vector<int> face;
                int g4loc = 0;
                if (is_hex) {
                    face.assign(HEX_FACE[side - 1], HEX_FACE[side - 1] + 4);
                    g4loc = HEX_G4[side - 1];
                } else if (is_tet) {
                    face.assign(TET_FACE[side - 1], TET_FACE[side - 1] + 3);
                    g4loc = TET_G4[side - 1];
                } else {
                    face.assign(PRI_FACE[side - 1], PRI_FACE[side - 1] + 4);
                    if (face.back() == 0) face.pop_back();
                    g4loc = PRI_G4[side - 1];
                }
                std::vector<double> p;
                for (size_t i = 0; i < face.size(); ++i)
                    p.push_back(i < l.RLOADi.size() ? l.RLOADi[i] : 0.0);

                // make sure the face order gives an outward normal
                const auto cc = corner_coords(model, e.NODIN, nc);
                if (!cc.empty()) {
                    const Vec3 p1 = cc[face[0] - 1], p2 = cc[face[1] - 1],
                               p3 = cc[face[2] - 1];
                    const Vec3 n = cross(sub(p2, p1), sub(p3, p1));
                    Vec3 cen{}, fcnt{};
                    for (const auto &c : cc) {
                        cen.x += c.x / (double)nc;
                        cen.y += c.y / (double)nc;
                        cen.z += c.z / (double)nc;
                    }
                    for (int f : face) {
                        fcnt.x += cc[f - 1].x / (double)face.size();
                        fcnt.y += cc[f - 1].y / (double)face.size();
                        fcnt.z += cc[f - 1].z / (double)face.size();
                    }
                    if (dot(n, sub(fcnt, cen)) < 0.0) {
                        if (face.size() == 4) {
                            std::swap(face[1], face[3]);
                            std::swap(p[1], p[3]);
                        } else {
                            std::swap(face[1], face[2]);
                            std::swap(p[1], p[2]);
                        }
                    }
                }
                const long g1 = ext_node(model, e.NODIN[face[0] - 1]);
                const long g4 = g4loc
                    ? ext_node(model, e.NODIN[g4loc - 1]) : 0;
                write_pload4(sid, eid, p, g1, g4);
                ++n_pload4;
            } else {
                warn_skipped(st, "BEUSLO on beam/other element");
            }
        }
    }

    // Sesam sets (GSETMEMB + TDSETNAM): kept as informational comments;
    // member ids are external numbers for node sets
    if (!model.sets.empty()) {
        out << "$\n$ Sesam sets (GSETMEMB/TDSETNAM, informational only):\n";
        for (const auto &kv : model.sets) {
            out << "$ Set " << kv.first;
            auto ni = model.set_names.find(kv.first);
            if (ni != model.set_names.end())
                out << " \"" << ni->second << "\"";
            out << (kv.second.first == 2 ? " (elements):" : " (nodes):");
            size_t n = 0;
            for (long m : kv.second.second) {
                if (n % 12 == 0) out << "\n$  ";
                out << " "
                    << (kv.second.first == 1 ? ext_node(model, m) : m);
                ++n;
            }
            out << "\n";
        }
    }

    out << "ENDDATA\n";
    out.close();

    // --------------------------------------------------------------
    // Report
    // --------------------------------------------------------------
    std::cout << "fem2bdf: converted '" << in_path << "' -> '" << out_path << "'\n"
              << "  cards read      : " << st.cards_read << "\n"
              << "  nodes (GRID)    : " << n_grid << "\n"
              << "  materials (MAT1): " << model.mats.size() + st.default_mats << "\n"
              << "  elements written: CQUAD4=" << n_quad
              << " CTRIA3=" << n_tria
              << " CBEAM=" << n_beam
              << " CROD=" << n_rod << "\n"
              << "                    CHEXA=" << n_hex
              << " CTETRA=" << n_tet
              << " CPENTA=" << n_pri
              << " CBUSH=" << n_bush << "\n"
              << "                    CONM1=" << n_conm1
              << " CMASS2=" << n_cmass2 << "\n"
              << "  load cases      : " << load_cases.size()
              << " (nodal: " << model.loads.size()
              << ", surface->PLOAD4: " << n_pload4 << ")\n"
              << "  SPC1 groups     : " << spc_groups.size()
              << (n_spc ? " (+" + std::to_string(n_spc) + " SPC entries)" : "")
              << "\n";

    if (!st.skipped_eltypes.empty()) {
        std::cout << "  skipped elements (unsupported type):\n";
        for (const auto &kv : st.skipped_eltypes)
            std::cout << "    type " << kv.first
                      << " (" << eltype_name(static_cast<fe::el_types>(kv.first))
                      << "): " << kv.second << "\n";
    }
    if (st.missing_geono)
        std::cout << "  WARNING: " << st.missing_geono
                  << " propert(ies) reference unsupported/missing GEONO cards\n";
    if (st.no_geono)
        std::cout << "  WARNING: " << st.no_geono
                  << " element(s) without geometry reference (GEONO=0) skipped\n";
    if (st.default_thk)
        std::cout << "  WARNING: " << st.default_thk
                  << " shell propert(ies) without GELTH thickness (or TH<=0);"
                     " zero-thickness placeholder PSHELL emitted - solver"
                     " will fail until corrected\n";
    if (st.no_prop)
        std::cout << "  WARNING: " << st.no_prop
                  << " element(s) skipped; their GEONO has no usable"
                     " section/property card\n";
    if (st.missing_mat)
        std::cout << "  WARNING: " << st.missing_mat
                  << " element(s) reference missing MISOSEL material\n";
    if (st.default_mats)
        std::cout << "  WARNING: " << st.default_mats
                  << " material(s) without MISOSEL; all-zero placeholder"
                     " MAT1 emitted - solver will fail until corrected\n";
    if (st.shcen_adjust)
        std::cout << "  NOTE: " << st.shcen_adjust
                  << " beam offset(s) shifted from centroid to shear centre"
                     " (PBEAML shapes are drawn from the shear centre)\n";
    if (st.gecc_local)
        std::cout << "  NOTE: " << st.gecc_local
                  << " GECC local eccentricit(ies) rotated to global\n";
    if (st.hinge_applied)
        std::cout << "  NOTE: " << st.hinge_applied
                  << " beam end(s) with pin flags from BELFIX\n";
    if (st.hinge_partial)
        std::cout << "  WARNING: " << st.hinge_partial
                  << " BELFIX dof(s) with fractional fixation/spring rate;"
                     " treated as fully connected (pin flags are binary)\n";
    if (st.hinge_spring)
        std::cout << "  WARNING: " << st.hinge_spring
                  << " BELFIX hinge(s) of OPT=SPRING; elastic rates dropped\n";
    if (st.hinge_local)
        std::cout << "  WARNING: " << st.hinge_local
                  << " BELFIX with TRANO (local system) ignored\n";
    if (st.missing_hinge)
        std::cout << "  WARNING: " << st.missing_hinge
                  << " beam end(s) reference unknown BELFIX FIXNO\n";
    if (st.morsmel_q)
        std::cout << "  WARNING: " << st.morsmel_q
                  << " MORSMEL material(s) with direction vector Q dropped\n";
    if (st.n_mpc)
        std::cout << "  NOTE: " << st.n_mpc
                  << " MPC card(s) written from BLDEP (set " << mpc_sid
                  << ")\n";
    if (st.missing_ldep)
        std::cout << "  WARNING: " << st.missing_ldep
                  << " dependent dof(s) (BNBCD=3) without matching BLDEP\n";
    if (st.missing_spring)
        std::cout << "  WARNING: " << st.missing_spring
                  << " GSPR element(s) without MGSPRNG property skipped\n";
    if (st.missing_mass)
        std::cout << "  WARNING: " << st.missing_mass
                  << " GMAS element(s) without MGMASS property skipped\n";
    if (st.spring_coupling)
        std::cout << "  WARNING: " << st.spring_coupling
                  << " spring(s) with off-diagonal stiffness; coupling dropped\n";
    if (st.spring_local)
        std::cout << "  WARNING: " << st.spring_local
                  << " spring/mass element(s) in local systems; "
                     "matrix used as global\n";
    if (st.accel_disp)
        std::cout << "  WARNING: " << st.accel_disp
                  << " prescribed-acceleration ent(ies) treated as fixed\n";
    if (st.missing_prescribed)
        std::cout << "  WARNING: " << st.missing_prescribed
                  << " prescribed dof(s) without BNDISPL; value 0 used\n";
    if (st.complx_disp)
        std::cout << "  WARNING: " << st.complx_disp
                  << " complex prescribed displacement(s); imaginary part dropped\n";
    if (st.missing_elem_load)
        std::cout << "  WARNING: " << st.missing_elem_load
                  << " surface load(s) reference unknown element(s)\n";
    if (st.comp_form_load)
        std::cout << "  WARNING: " << st.comp_form_load
                  << " component-form surface load(s) (LOTYP=2) skipped\n";
    if (st.hydro_dummy)
        std::cout << "  WARNING: " << st.hydro_dummy
                  << " hydro-pressure marker(s) (LOTYP=3) skipped\n";
    if (st.follower_load)
        std::cout << "  WARNING: " << st.follower_load
                  << " non-conservative surface load(s) converted as normal pressure\n";
    if (st.complx_load)
        std::cout << "  WARNING: " << st.complx_load
                  << " complex surface load(s); imaginary part dropped\n";
    if (st.side_mid)
        std::cout << "  WARNING: " << st.side_mid
                  << " mid-surface pressure(s) (SIDE=2) treated as top face\n";
    if (st.bad_side)
        std::cout << "  WARNING: " << st.bad_side
                  << " surface load(s) with unusable SIDE flag skipped\n";
    if (st.trapezoid_bar)
        std::cout << "  WARNING: " << st.trapezoid_bar
                  << " trapezoidal GBARM section(s) converted to equivalent PBEAM\n";
    if (st.usyi_offset)
        std::cout << "  WARNING: " << st.usyi_offset
                  << " GUSYI section(s) with off-centre web; offsets dropped\n";
    if (st.curvature_ignored)
        std::cout << "  WARNING: " << st.curvature_ignored
                  << " section(s) with inside curvature; fillet radii dropped\n";
    if (st.tonp_as_i)
        std::cout << "  NOTE: " << st.tonp_as_i
                  << " GTONP section(s) converted to I (plate as bottom flange)\n";
    if (!st.skipped_cards.empty()) {
        std::cout << "  skipped cards:\n";
        for (const auto &kv : st.skipped_cards)
            std::cout << "    " << kv.first << ": " << kv.second << "\n";
    }
    return 0;
}
