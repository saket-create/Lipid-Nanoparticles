"""
LNP PROCESS SIMULATION — COMPOSITION & MASS MODULE

Reads Section B (feed streams: aqueous, organic, dilution buffer, cryo feed)
and Section C (final product stream) from lnp_input.txt.

Core equation:  mass_flow (mg/min) = concentration (mg/mL) x flow_rate (mL/min)

CALCULATION MODES (detected automatically):
  FORWARD    Section B filled, Section C blank.
             Feed concentrations propagate through all process steps.
             Derived streams: mixer output, TFF retentate, TFF permeate, final product.

  BACKWARD   Section C filled, Section B blank.
             Final product composition is back-calculated to feed streams using
             TFF retention fractions in reverse and flow rates from the flow solver.

  BOTH       Sections B and C both filled.
             Forward calculation is performed from feed data.
             Final product values are then compared against Section C inputs.
             Any discrepancy triggers a WARNING.

TFF separation uses realistic efficiency fractions (not ideal 100% cuts).
The flow rate model (Section A / lnp_simulation.py) is always the authority.
Composition data never modifies flow rates.
"""

import os
from datetime import datetime

AQUEOUS_COMPS  = ["mRNA", "Citric_Acid", "Sodium_Citrate", "EDTA", "Water"]
ORGANIC_COMPS  = ["Ionisable_Lipid", "DSPC", "Cholesterol", "PEG_Lipid", "Ethanol"]
DIL_BUF_COMPS  = ["NaCl", "Na2HPO4", "NaH2PO4", "Water"]
CRYO_COMPS     = ["Trehalose", "Water"]

FP_COMPS = sorted(
    set(AQUEOUS_COMPS) | {"Ionisable_Lipid","DSPC","Cholesterol","PEG_Lipid","Ethanol"}
    | set(CRYO_COMPS)
)

ALL_COMP_NAMES = set(
    AQUEOUS_COMPS + ORGANIC_COMPS + DIL_BUF_COMPS + CRYO_COMPS
)

PRIMARY_STREAMS = ["aqueous", "organic", "dilution_buffer", "cryo_feed"]
STREAM_COMPS = {
    "aqueous"        : AQUEOUS_COMPS,
    "organic"        : ORGANIC_COMPS,
    "dilution_buffer": DIL_BUF_COMPS,
    "cryo_feed"      : CRYO_COMPS,
}

TFF_RETENTION_FRACTION = {
    "mRNA"            : 0.95,
    "Ionisable_Lipid" : 0.98,
    "DSPC"            : 0.98,
    "Cholesterol"     : 0.98,
    "PEG_Lipid"       : 0.98,
    "Ethanol"         : 0.05,
    "Citric_Acid"     : 0.02,
    "Sodium_Citrate"  : 0.02,
    "EDTA"            : 0.02,
    "NaCl"            : 0.02,
    "Na2HPO4"         : 0.02,
    "NaH2PO4"         : 0.02,
    "Water"           : 0.10,
    "Trehalose"       : 0.00,
}

COMP_DISPLAY = {
    "mRNA"            : "mRNA",
    "Citric_Acid"     : "Citric Acid",
    "Sodium_Citrate"  : "Sodium Citrate",
    "EDTA"            : "EDTA",
    "Water"           : "Water",
    "Ionisable_Lipid" : "Ionisable Lipid (SM-102/MC3)",
    "DSPC"            : "DSPC (Phospholipid)",
    "Cholesterol"     : "Cholesterol",
    "PEG_Lipid"       : "PEG-DMG / PEG-Lipid",
    "Ethanol"         : "Ethanol",
    "NaCl"            : "Sodium Chloride (NaCl)",
    "Na2HPO4"         : "Dibasic Sodium Phosphate",
    "NaH2PO4"         : "Monobasic Sodium Phosphate",
    "Trehalose"       : "Trehalose / Sucrose (Cryoprotectant)",
}

STREAM_TITLES = {
    "aqueous"        : "AQUEOUS PHASE  (mRNA + Buffer, Pump P-101)",
    "organic"        : "ORGANIC PHASE  (Lipids in Ethanol, Pump P-102)",
    "dilution_buffer": "DILUTION BUFFER  (PBS / Buffer, T-105)",
    "cryo_feed"      : "CRYOPROTECTANT FEED  (Trehalose/Sucrose, T-100)",
    "mixer_out"      : "MIXER OUTPUT  (Crude LNP Dispersion, M-101)",
    "tff_retentate"  : "TFF RETENTATE  (Purified Concentrated LNP, UF-101)",
    "tff_permeate"   : "TFF PERMEATE  (Waste Stream, UF-101)",
    "final_product"  : "FINAL PRODUCT  (Sterile LNP Formulation, ST-101)",
}

STREAM_FLOW_VAR = {
    "aqueous"        : "flow_aqueous_mL_min",
    "organic"        : "flow_organic_mL_min",
    "dilution_buffer": "dilution_flow_mL_min",
    "cryo_feed"      : "cryo_flow_mL_min",
    "mixer_out"      : "mixer_flow_mL_min",
    "tff_retentate"  : "retentate_mL_min",
    "tff_permeate"   : "permeate_mL_min",
    "final_product"  : "final_product_mL_min",
}

STREAM_ORDER = [
    "aqueous","organic","dilution_buffer","cryo_feed",
    "mixer_out","tff_retentate","tff_permeate","final_product",
]

MASS_BAL_TOL  = 1e-3
THREE_WAY_TOL = 1e-3
CONC_LIMIT    = 1100.0
VERIFY_TOL    = 0.05  # 5% relative tolerance for back-calculation verification


def _add(a, b):
    if a is None and b is None:
        return None
    return (a or 0.0) + (b or 0.0)


def _build_flows(sol, params):
    flows = {k: v for k, v in sol.items() if isinstance(v, float)}
    flows["cryo_flow_mL_min"] = params.get("cryo_flow_mL_min", 0.0)
    return flows


def fv(v):
    return f"{v:.4f}" if v is not None else "—"


def _r(comp):
    return TFF_RETENTION_FRACTION.get(comp, 0.02)


# ── Validation helpers ─────────────────────────────────────────────────────────

def _check_negatives(stream, comp_data, issues):
    for comp, r in comp_data.items():
        disp = COMP_DISPLAY.get(comp, comp)
        for field, val in [("concentration", r.get("conc")), ("mass flow", r.get("mass"))]:
            if val is not None and val < 0:
                issues.append({"level":"ERROR",
                    "msg": f"Negative {field} for {disp} in '{stream}' "
                           f"(value = {val:.4f}). All values must be positive."})

def _check_total_conc(stream, comp_data, issues):
    vals = [r.get("conc") for r in comp_data.values() if r.get("conc") is not None]
    if not vals:
        return
    total = sum(vals)
    if total > CONC_LIMIT:
        issues.append({"level":"WARNING",
            "msg": f"Total concentration in '{stream}' = {total:.1f} mg/mL exceeds "
                   f"{CONC_LIMIT:.0f} mg/mL. Pure water is ~1000 mg/mL. "
                   f"The solute concentrations in this stream appear too high. "
                   f"Please check your input values."})

def _check_three_way(stream, comp, conc, mass, flow_rate, issues):
    implied = conc * flow_rate
    diff = abs(mass - implied)
    if diff > THREE_WAY_TOL:
        disp = COMP_DISPLAY.get(comp, comp)
        issues.append({"level":"WARNING",
            "msg": f"Inconsistency for {disp} in '{stream}': "
                   f"given mass = {mass:.4f} mg/min, but "
                   f"conc ({conc:.4f} mg/mL) x flow ({flow_rate:.4f} mL/min) "
                   f"= {implied:.4f} mg/min. "
                   f"The flow rate model is the authority. "
                   f"Using {implied:.4f} mg/min for all calculations."})
    return implied


# ── A: Forward path — solve primary feed streams ───────────────────────────────

def solve_feed_streams(comp_known, flows, issues):
    """
    Solves the four primary (user-input) feed streams.
    For each component: conc given → mass = conc x flow; mass given → conc = mass / flow.
    If both given, three-way check is applied. Flow rate is always the authority.
    Returns dict: stream → {comp: {conc, mass, src}}
    """
    results = {}
    for stream in PRIMARY_STREAMS:
        comps     = STREAM_COMPS[stream]
        flow_var  = STREAM_FLOW_VAR[stream]
        flow_rate = flows.get(flow_var)
        stream_result = {}

        for comp in comps:
            conc = comp_known.get((stream, comp, "conc_mg_mL"))
            mass = comp_known.get((stream, comp, "mass_mg_min"))
            disp = COMP_DISPLAY.get(comp, comp)

            if conc is not None and mass is not None and flow_rate:
                implied = _check_three_way(stream, comp, conc, mass, flow_rate, issues)
                stream_result[comp] = {"conc": conc, "mass": implied,
                    "src": "conc given" if abs(mass-implied) <= THREE_WAY_TOL
                           else "conc given (mass overridden by flow model)"}

            elif conc is not None and flow_rate:
                stream_result[comp] = {"conc": conc, "mass": conc*flow_rate, "src": "conc given"}

            elif mass is not None and flow_rate:
                stream_result[comp] = {"conc": mass/flow_rate, "mass": mass, "src": "mass given"}

            elif conc is not None:
                stream_result[comp] = {"conc": conc, "mass": None, "src": "conc given, no flow"}

            elif mass is not None:
                stream_result[comp] = {"conc": None, "mass": mass, "src": "mass given, no flow"}

            else:
                stream_result[comp] = {"conc": None, "mass": None, "src": "not provided"}

        _check_negatives(stream, stream_result, issues)
        _check_total_conc(stream, stream_result, issues)
        results[stream] = stream_result

    return results


# ── B: Forward path — propagate derived streams ────────────────────────────────

def propagate_forward(feed, flows):
    """
    Propagates composition from feed streams through the process.
    Uses mass conservation and TFF_RETENTION_FRACTION for TFF step.
    Returns dict of four derived streams: mixer_out, tff_retentate, tff_permeate, final_product
    """
    aq  = feed["aqueous"]
    org = feed["organic"]
    dil = feed["dilution_buffer"]
    cry = feed["cryo_feed"]

    mixer_flow = flows.get("mixer_flow_mL_min")
    mixer_comps = set(AQUEOUS_COMPS) | set(ORGANIC_COMPS)
    mixer_out = {}
    for comp in sorted(mixer_comps):
        m = _add(aq.get(comp,{}).get("mass"), org.get(comp,{}).get("mass"))
        c = m/mixer_flow if (m is not None and mixer_flow) else None
        mixer_out[comp] = {"conc":c, "mass":m, "src":"calculated"}

    pre_tff = {}
    for comp in sorted(mixer_comps):
        pre_tff[comp] = _add(pre_tff.get(comp), mixer_out[comp]["mass"])
    for comp in DIL_BUF_COMPS:
        pre_tff[comp] = _add(pre_tff.get(comp), dil.get(comp,{}).get("mass"))

    ret_flow  = flows.get("retentate_mL_min")
    perm_flow = flows.get("permeate_mL_min")
    all_pre_comps = sorted(set(mixer_comps) | set(DIL_BUF_COMPS))
    tff_ret  = {}
    tff_perm = {}
    for comp in all_pre_comps:
        m_total = pre_tff.get(comp)
        if m_total is None:
            tff_ret[comp]  = {"conc":None,"mass":None,"src":"calculated"}
            tff_perm[comp] = {"conc":None,"mass":None,"src":"calculated"}
            continue
        frac   = _r(comp)
        m_ret  = m_total * frac
        m_perm = m_total * (1.0 - frac)
        tff_ret[comp]  = {"conc": m_ret /ret_flow  if ret_flow  else None, "mass":m_ret,  "src":"calculated"}
        tff_perm[comp] = {"conc": m_perm/perm_flow if perm_flow else None, "mass":m_perm, "src":"calculated"}

    fp_flow = flows.get("final_product_mL_min")
    fp_comps = sorted(set(tff_ret.keys()) | set(CRYO_COMPS))
    fp = {}
    for comp in fp_comps:
        m = _add(tff_ret.get(comp,{}).get("mass"), cry.get(comp,{}).get("mass"))
        c = m/fp_flow if (m is not None and fp_flow) else None
        fp[comp] = {"conc":c,"mass":m,"src":"calculated"}

    return {"mixer_out":mixer_out,"tff_retentate":tff_ret,"tff_permeate":tff_perm,"final_product":fp}


# ── C: Backward path — back-calculate from final product ──────────────────────

def back_calculate(fp_known, flows, issues):
    """
    Given final product concentrations/masses (Section C), back-calculates
    the composition of upstream feed streams using TFF retention fractions in reverse.

    Back-calculation rules (component-specific):
      mRNA, Citric_Acid, Sodium_Citrate, EDTA  ->  exclusive to aqueous phase.
        aqueous_mass = fp_mass / TFF_retention_fraction
        aqueous_conc = aqueous_mass / aqueous_flow_rate

      Ionisable_Lipid, DSPC, Cholesterol, PEG_Lipid, Ethanol  ->  exclusive to organic.
        organic_mass = fp_mass / TFF_retention_fraction
        organic_conc = organic_mass / organic_flow_rate

      Trehalose  ->  added in cryo feed, not in TFF feed.
        cryo_mass = fp_mass  (all trehalose in product came from cryo feed)
        cryo_conc = cryo_mass / cryo_flow_rate

      Water  ->  present in aqueous, organic, dil buffer, and cryo feed.
        Cannot be unambiguously back-calculated without more data.
        Reported as indeterminate.

    Returns dict of all streams.
    """
    fp_flow   = flows.get("final_product_mL_min")
    ret_flow  = flows.get("retentate_mL_min")
    perm_flow = flows.get("permeate_mL_min")
    aq_flow   = flows.get("flow_aqueous_mL_min")
    org_flow  = flows.get("flow_organic_mL_min")
    cryo_flow = flows.get("cryo_flow_mL_min")

    AQUEOUS_ONLY = {"mRNA", "Citric_Acid", "Sodium_Citrate", "EDTA"}
    ORGANIC_ONLY = {"Ionisable_Lipid", "DSPC", "Cholesterol", "PEG_Lipid", "Ethanol"}
    CRYO_ONLY    = {"Trehalose"}

    fp_comps = set(AQUEOUS_COMPS) | set(ORGANIC_COMPS) | set(CRYO_COMPS)

    # Step 1: Solve final product stream
    fp = {}
    for comp in sorted(fp_comps):
        conc = fp_known.get(("final_product", comp, "conc_mg_mL"))
        mass = fp_known.get(("final_product", comp, "mass_mg_min"))
        if conc is not None and mass is not None and fp_flow:
            implied = _check_three_way("final_product", comp, conc, mass, fp_flow, issues)
            fp[comp] = {"conc": conc, "mass": implied, "src": "user input"}
        elif conc is not None and fp_flow:
            fp[comp] = {"conc": conc, "mass": conc * fp_flow, "src": "user input"}
        elif mass is not None and fp_flow:
            fp[comp] = {"conc": mass / fp_flow, "mass": mass, "src": "user input"}
        elif conc is not None:
            fp[comp] = {"conc": conc, "mass": None, "src": "user input"}
        elif mass is not None:
            fp[comp] = {"conc": None, "mass": mass, "src": "user input"}
        else:
            fp[comp] = {"conc": None, "mass": None, "src": "not provided"}

    _check_negatives("final_product", fp, issues)
    _check_total_conc("final_product", fp, issues)

    # Step 2: Back-calculate each stream
    aq  = {c: {"conc": None, "mass": None, "src": "back-calculated"} for c in AQUEOUS_COMPS}
    org = {c: {"conc": None, "mass": None, "src": "back-calculated"} for c in ORGANIC_COMPS}
    cry = {c: {"conc": None, "mass": None, "src": "back-calculated"} for c in CRYO_COMPS}
    tff_ret  = {}
    tff_perm = {}

    for comp in sorted(fp_comps):
        fp_mass = fp.get(comp, {}).get("mass")
        if fp_mass is None:
            continue

        frac = _r(comp)

        if comp in CRYO_ONLY:
            # Trehalose entirely from cryo feed
            cry[comp] = {"conc": fp_mass / cryo_flow if cryo_flow else None,
                         "mass": fp_mass, "src": "back-calculated"}
            tff_ret[comp] = {"conc": None, "mass": None, "src": "not in TFF feed"}
            tff_perm[comp] = {"conc": None, "mass": None, "src": "not in TFF feed"}
            continue

        if comp == "Water":
            # Water comes from multiple sources — indeterminate
            aq["Water"]  = {"conc": None, "mass": None, "src": "indeterminate (multiple sources)"}
            cry["Water"] = {"conc": None, "mass": None, "src": "indeterminate (multiple sources)"}
            ret_m = fp_mass
            perm_m = ret_m / frac * (1.0 - frac) if frac > 0 else None
            tff_ret["Water"]  = {"conc": ret_m / ret_flow if ret_flow else None,
                                  "mass": ret_m, "src": "back-calculated (partial)"}
            tff_perm["Water"] = {"conc": perm_m / perm_flow if (perm_m and perm_flow) else None,
                                  "mass": perm_m, "src": "back-calculated (partial)"}
            continue

        # For exclusive components: fp_mass = retentate_mass (cryo doesn't affect lipids/mRNA)
        if frac == 0:
            continue
        pre_tff_mass = fp_mass / frac
        perm_mass    = pre_tff_mass * (1.0 - frac)

        tff_ret[comp]  = {"conc": fp_mass / ret_flow if ret_flow else None,
                           "mass": fp_mass, "src": "back-calculated"}
        tff_perm[comp] = {"conc": perm_mass / perm_flow if perm_flow else None,
                           "mass": perm_mass, "src": "back-calculated"}

        if comp in AQUEOUS_ONLY and aq_flow:
            aq[comp] = {"conc": pre_tff_mass / aq_flow, "mass": pre_tff_mass,
                        "src": "back-calculated"}

        elif comp in ORGANIC_ONLY and org_flow:
            org[comp] = {"conc": pre_tff_mass / org_flow, "mass": pre_tff_mass,
                         "src": "back-calculated"}

    # Dilution buffer cannot be back-calculated from final product alone
    dil = {c: {"conc": None, "mass": None,
               "src": "indeterminate (dil buffer not separable from final product)"}
           for c in DIL_BUF_COMPS}

    # Mixer output
    mixer_flow = flows.get("mixer_flow_mL_min")
    mixer_out = {}
    for comp in sorted(set(AQUEOUS_COMPS) | set(ORGANIC_COMPS)):
        m = _add(aq.get(comp, {}).get("mass"), org.get(comp, {}).get("mass"))
        c = m / mixer_flow if (m is not None and mixer_flow) else None
        mixer_out[comp] = {"conc": c, "mass": m, "src": "back-calculated"}

    _check_negatives("aqueous",       aq,      issues)
    _check_negatives("organic",       org,     issues)
    _check_negatives("tff_retentate", tff_ret, issues)

    return {
        "aqueous": aq, "organic": org, "dilution_buffer": dil, "cryo_feed": cry,
        "mixer_out": mixer_out, "tff_retentate": tff_ret,
        "tff_permeate": tff_perm, "final_product": fp,
    }


def verify_final_product(fp_calculated, fp_user_known, fp_flow, issues):
    """
    When both feed (Section B) and final product (Section C) data are given,
    compare the forward-calculated final product against the user's Section C values.
    Any relative difference > VERIFY_TOL (5%) is flagged as a WARNING.
    Returns list of verification result dicts.
    """
    results = []
    for comp in sorted(set(fp_calculated) | set(k[1] for k in fp_user_known if k[0]=="final_product")):
        for field, label in [("conc_mg_mL","conc"),("mass_mg_min","mass")]:
            user_val = fp_user_known.get(("final_product", comp, field))
            if user_val is None:
                continue
            calc_val = fp_calculated.get(comp, {}).get(label)
            if calc_val is None:
                continue
            diff  = abs(calc_val - user_val)
            rel   = diff / abs(user_val) if user_val != 0 else float("inf")
            disp  = COMP_DISPLAY.get(comp, comp)
            unit  = "mg/mL" if label == "conc" else "mg/min"
            passed = rel <= VERIFY_TOL
            if not passed:
                issues.append({"level":"WARNING",
                    "msg": f"Final product mismatch for {disp} ({label}): "
                           f"calculated = {calc_val:.4f} {unit}, "
                           f"user given = {user_val:.4f} {unit}, "
                           f"relative difference = {rel*100:.1f}% "
                           f"(tolerance = {VERIFY_TOL*100:.0f}%)."})
            results.append({
                "comp": comp, "field": label, "disp": disp,
                "calculated": calc_val, "user_given": user_val,
                "diff": diff, "rel_pct": rel*100, "passed": passed, "unit": unit,
            })
    return results


# ── E: Component-wise mass balance ─────────────────────────────────────────────

def check_mass_balance(feed, tff_ret, tff_perm):
    """
    At TFF level: aq + org + dil = retentate + permeate for every component.
    """
    results = []
    aq, org, dil = feed["aqueous"], feed["organic"], feed["dilution_buffer"]
    all_comps = sorted(
        set(aq)|set(org)|set(dil)|set(tff_ret)|set(tff_perm)
    )
    for comp in all_comps:
        m_in  = sum((s.get(comp,{}).get("mass") or 0.0) for s in [aq,org,dil])
        m_out = sum((s.get(comp,{}).get("mass") or 0.0) for s in [tff_ret,tff_perm])
        has_in  = any(s.get(comp,{}).get("mass") is not None for s in [aq,org,dil])
        has_out = any(s.get(comp,{}).get("mass") is not None for s in [tff_ret,tff_perm])
        if not (has_in and has_out):
            continue
        err    = abs(m_in - m_out)
        pct    = err/m_in*100 if m_in > 0 else 0.0
        results.append({
            "comp": comp, "disp": COMP_DISPLAY.get(comp, comp),
            "mass_in": m_in, "mass_out": m_out,
            "error": err, "pct": pct, "passed": err < MASS_BAL_TOL,
        })
    return results


# ── F: Report writer ───────────────────────────────────────────────────────────

def write_report(folder, run_id, all_results, flows, issues, mb_results,
                 fp_verify, calc_mode):
    now = datetime.now().strftime("%Y-%m-%d  %H:%M:%S")
    L   = []
    def a(s=""): L.append(s)

    a("  LNP COMPOSITION & MASS REPORT")
    a(f"  Run ID : {run_id:03d}    |    Timestamp : {now}")
    a(f"  Calculation mode : {calc_mode}")
    a()
    a("  mass_flow (mg/min) = concentration (mg/mL) x flow_rate (mL/min)")
    a("  TFF separation uses realistic efficiency fractions (see table below).")
    a()

    errors   = [i for i in issues if i["level"] == "ERROR"]
    warnings = [i for i in issues if i["level"] == "WARNING"]

    if errors:
        a("  ERRORS")
        for e in errors: a(f"    {e['msg']}")
        a()
    if warnings:
        a("  WARNINGS")
        for w in warnings: a(f"    {w['msg']}")
        a()

    a("  TFF SEPARATION EFFICIENCIES")
    a(f"  {'Component':<35} {'Retained (retentate)':>22} {'Removed (permeate)':>20}")
    for comp in sorted(TFF_RETENTION_FRACTION):
        if comp == "Trehalose": continue
        frac = TFF_RETENTION_FRACTION[comp]
        a(f"  {COMP_DISPLAY.get(comp,comp):<35} {frac*100:>21.0f}% {(1-frac)*100:>19.0f}%")
    a()

    for stream in STREAM_ORDER:
        comps  = all_results.get(stream, {})
        fvar   = STREAM_FLOW_VAR[stream]
        fr     = flows.get(fvar)
        fr_str = f"{fr:.4f} mL/min" if fr is not None else "not available"
        a(f"  {STREAM_TITLES[stream]}")
        a(f"  Flow rate : {fr_str}")
        a()
        a(f"  {'Component':<35} {'Conc (mg/mL)':>14} {'Mass Flow (mg/min)':>20}  Source")
        if not comps:
            a("  (no data)")
        else:
            for comp, r in comps.items():
                src = r.get("src","")
                tag = "(user input)" if "user" in src or "given" in src else "(back-calc)" if "back" in src else "(calculated)"
                a(f"  {COMP_DISPLAY.get(comp,comp):<35} {fv(r.get('conc')):>14} {fv(r.get('mass')):>20}  {tag}")
        a()

    if fp_verify:
        a("  FINAL PRODUCT VERIFICATION  (Section B vs Section C)")
        a("  Calculated from feed data vs your Section C input values")
        a(f"  Tolerance : {VERIFY_TOL*100:.0f}% relative difference")
        a()
        a(f"  {'Component':<35} {'Field':<12} {'Calculated':>12} {'User Given':>12} {'Rel Diff':>10}  Status")
        for r in fp_verify:
            status = "PASS" if r["passed"] else "FAIL"
            a(f"  {r['disp']:<35} {r['field']:<12} {r['calculated']:>12.4f} "
              f"{r['user_given']:>12.4f} {r['rel_pct']:>9.1f}%  {status}")
        n_fail = sum(1 for r in fp_verify if not r["passed"])
        a()
        if n_fail == 0:
            a("  RESULT: Final product values are CONSISTENT with feed data.")
        else:
            a(f"  RESULT: {n_fail} value(s) INCONSISTENT. Review feed or product data.")
        a()

    if mb_results:
        a("  COMPONENT-WISE MASS BALANCE  (TFF level: inputs = retentate + permeate)")
        a()
        a(f"  {'Component':<35} {'Mass In':>12} {'Mass Out':>12} {'Error':>10}  Status")
        all_pass = True
        for r in mb_results:
            status = "PASS" if r["passed"] else f"FAIL ({r['error']:.4f} mg/min, {r['pct']:.2f}%)"
            if not r["passed"]: all_pass = False
            a(f"  {r['disp']:<35} {r['mass_in']:>12.4f} {r['mass_out']:>12.4f} {r['error']:>10.4f}  {status}")
        a()
        a(f"  Mass balance: {'ALL COMPONENTS BALANCED' if all_pass else 'SOME COMPONENTS FAILED'}")
        a()

    a("  END OF COMPOSITION REPORT")
    out_path = os.path.join(folder, "composition_report.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    return out_path


# ── G: Main solver ─────────────────────────────────────────────────────────────

def run_composition_from_solver(comp_known, sol, params, folder, run_id):
    """
    Main entry point called by lnp_simulation.py.
    comp_known contains both Section B (feed) and Section C (final product) entries.
    """
    print("  LNP COMPOSITION & MASS MODULE")
    print()

    flows  = _build_flows(sol, params)
    issues = []

    has_feed_b = any(k[0] in PRIMARY_STREAMS for k in comp_known)
    has_fp_c   = any(k[0] == "final_product" for k in comp_known)

    if has_feed_b and has_fp_c:
        calc_mode = "BOTH (forward from feed, then verify against final product)"
    elif has_feed_b:
        calc_mode = "FORWARD (feed streams → derived streams)"
    elif has_fp_c:
        calc_mode = "BACKWARD (final product → back-calculate feed streams)"
    else:
        print("  No composition data found in Section B or C of the input file.")
        return

    print(f"  Mode: {calc_mode}")
    print()

    fp_verify  = []
    mb_results = []

    if has_feed_b:
        feed    = solve_feed_streams(comp_known, flows, issues)
        derived = propagate_forward(feed, flows)
        all_results = {**feed, **derived}

        if has_fp_c:
            fp_verify = verify_final_product(
                derived["final_product"], comp_known,
                flows.get("final_product_mL_min"), issues
            )

        mb_results = check_mass_balance(
            feed, derived["tff_retentate"], derived["tff_permeate"]
        )

    else:
        all_results = back_calculate(comp_known, flows, issues)

    errs  = [i for i in issues if i["level"] == "ERROR"]
    warns = [i for i in issues if i["level"] == "WARNING"]

    if errs:
        print("  ERRORS:")
        for e in errs: print(f"    {e['msg']}")
        print()
    if warns:
        print("  WARNINGS:")
        for w in warns: print(f"    {w['msg']}")
        print()

    for stream in STREAM_ORDER:
        comps  = all_results.get(stream, {})
        fvar   = STREAM_FLOW_VAR[stream]
        fr     = flows.get(fvar)
        fr_str = f"{fr:.4f} mL/min" if fr is not None else "no flow rate"
        print(f"  {STREAM_TITLES[stream]}  [{fr_str}]")
        for comp, r in comps.items():
            c_str = f"{r['conc']:.4f} mg/mL"  if r.get("conc") is not None else "—"
            m_str = f"{r['mass']:.4f} mg/min" if r.get("mass") is not None else "—"
            print(f"    {COMP_DISPLAY.get(comp,comp):<38}  {c_str:<22}  {m_str}")
        print()

    if fp_verify:
        n_fail = sum(1 for r in fp_verify if not r["passed"])
        print(f"  Final Product Verification: {'ALL CONSISTENT' if n_fail==0 else f'{n_fail} MISMATCH(ES) FOUND'}")
        print()

    if mb_results:
        n_fail = sum(1 for r in mb_results if not r["passed"])
        print(f"  Component Mass Balance: {'ALL PASS' if n_fail==0 else f'{n_fail} COMPONENT(S) FAILED'}")
        print()

    out = write_report(folder, run_id, all_results, flows,
                       issues, mb_results, fp_verify, calc_mode)
    print(f"  Composition report saved: {out}")
    print()


# ── H: Standalone reader helpers ───────────────────────────────────────────────

def _read_comp_known_from_input(filepath="lnp_input.txt"):
    """Reads Section B and Section C composition lines from lnp_input.txt."""
    comp_known    = {}
    valid_streams = set(PRIMARY_STREAMS) | {"final_product"}
    valid_fields  = {"conc_mg_mL", "mass_mg_min"}
    valid_comps   = ALL_COMP_NAMES

    if not os.path.isfile(filepath):
        return comp_known
    with open(filepath, encoding="utf-8") as f:
        lines = f.readlines()
    for line in lines:
        parse = line.split("|")[0].strip()
        if not parse or "=" not in parse or "." not in parse:
            continue
        key, _, val = parse.partition("=")
        key, val = key.strip(), val.strip()
        if not val:
            continue
        parts = key.split(".")
        if len(parts) == 3:
            stream, comp, field = parts
            if stream in valid_streams and comp in valid_comps and field in valid_fields:
                try:
                    comp_known[(stream, comp, field)] = float(val)
                except ValueError:
                    pass
    return comp_known


def _read_flows_from_latest_report():
    base = "simulation_runs"
    if not os.path.isdir(base):
        return {}, None, 0
    runs = sorted([d for d in os.listdir(base) if d.startswith("run_")])
    if not runs:
        return {}, None, 0
    folder = os.path.join(base, runs[-1])
    try: run_id = int(runs[-1].split("_")[-1])
    except ValueError: run_id = 0
    report = os.path.join(folder, "output_report.txt")
    if not os.path.isfile(report):
        return {}, folder, run_id
    label_to_var = {
        "Aqueous Phase Flow Rate (into Mixer)"  : "flow_aqueous_mL_min",
        "Organic Phase Flow Rate (into Mixer)"  : "flow_organic_mL_min",
        "Mixer Output Flow Rate"                : "mixer_flow_mL_min",
        "Dilution Buffer Flow Rate"             : "dilution_flow_mL_min",
        "TFF Retentate Flow Rate (Purified LNP)": "retentate_mL_min",
        "TFF Permeate Flow Rate (Waste)"        : "permeate_mL_min",
        "Flow After Cryoprotectant Addition"    : "after_cryo_flow_mL_min",
        "Final Product Flow Rate"               : "final_product_mL_min",
        "Cryoprotectant Addition Flow Rate"     : "cryo_flow_mL_min",
        "Cryoprotectant Flow Rate"              : "cryo_flow_mL_min",
    }
    flows = {}
    with open(report, encoding="utf-8") as f:
        for line in f:
            for label, var in label_to_var.items():
                if label in line and var not in flows:
                    for token in line.strip().split():
                        try: flows[var] = float(token); break
                        except ValueError: continue
    return flows, folder, run_id


if __name__ == "__main__":
    print()
    print("  LNP COMPOSITION & MASS MODULE  (standalone)")
    print("  Reading from: lnp_input.txt")
    print("  Flow rates from: latest simulation_runs/run_XXX/output_report.txt")
    print()
    comp_known = _read_comp_known_from_input("lnp_input.txt")
    if not comp_known:
        print("  No composition data found in Section B or C of lnp_input.txt.")
        print("  Fill in concentrations or mass flows and re-run.")
    else:
        flows, folder, run_id = _read_flows_from_latest_report()
        if not flows:
            print("  No flow rate results found. Run lnp_simulation.py first.")
        else:
            run_composition_from_solver(
                comp_known, flows,
                {"cryo_flow_mL_min": flows.get("cryo_flow_mL_min", 0.0)},
                folder, run_id
            )
