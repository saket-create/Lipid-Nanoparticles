################################################################################
#  LNP MANUFACTURING PROCESS SIMULATION
#  Main simulation script.
#  Input file  : lnp_input.txt
#  User manual : lnp_user_manual.txt
#  Usage       : python lnp_simulation.py
################################################################################

import os, shutil
import numpy as np
from datetime import datetime

INPUT_FILE = "lnp_input.txt"

# ── Variable and parameter definitions ─────────────────────────────────────────

FLOW_VARIABLES = [
    "flow_aqueous_mL_min",
    "flow_organic_mL_min",
    "mixer_flow_mL_min",
    "dilution_flow_mL_min",
    "tff_feed_mL_min",
    "retentate_mL_min",
    "permeate_mL_min",
    "after_cryo_flow_mL_min",
    "final_product_mL_min",
]
ALL_VARIABLES  = FLOW_VARIABLES + ["FRR"]
PROCESS_PARAMS = [
    "dilution_ratio",
    "tff_vrr",
    "cryo_flow_mL_min",
    "sterile_filter_dead_vol_mL_min",
]

DISPLAY = {
    "flow_aqueous_mL_min"   : "Aqueous Phase Flow Rate (into Mixer)",
    "flow_organic_mL_min"   : "Organic Phase Flow Rate (into Mixer)",
    "FRR"                   : "Flow Rate Ratio, FRR (Aqueous / Organic)",
    "mixer_flow_mL_min"     : "Mixer Output Flow Rate",
    "dilution_flow_mL_min"  : "Dilution Buffer Flow Rate",
    "tff_feed_mL_min"       : "TFF Feed Flow Rate",
    "retentate_mL_min"      : "TFF Retentate Flow Rate (Purified LNP)",
    "permeate_mL_min"       : "TFF Permeate Flow Rate (Waste)",
    "after_cryo_flow_mL_min": "Flow After Cryoprotectant Addition",
    "final_product_mL_min"  : "Final Product Flow Rate",
}

# Knowability score: how likely is a user to know this value independently?
# Higher = more directly measured/set. Lower = more derived/indirect.
# Used to choose which variables to "assume unknown" in over-determined cases.
KNOWABILITY = {
    "flow_aqueous_mL_min"   : 0.95,  # Directly set by pump P-101
    "flow_organic_mL_min"   : 0.95,  # Directly set by pump P-102
    "FRR"                   : 0.85,  # Ratio of pump settings
    "mixer_flow_mL_min"     : 0.70,  # Sum of pump flows
    "dilution_flow_mL_min"  : 0.90,  # Set by dilution pump
    "tff_feed_mL_min"       : 0.60,  # Derived; rarely directly measured
    "retentate_mL_min"      : 0.50,  # TFF-dependent; variable
    "permeate_mL_min"       : 0.45,  # Derived from TFF; least directly known
    "after_cryo_flow_mL_min": 0.55,  # Intermediate derived variable
    "final_product_mL_min"  : 0.80,  # Target or measured output
}

# Realistic operating ranges for each variable [min, max]
# Used by the physical validation check.
REALISTIC_RANGES = {
    "flow_aqueous_mL_min"   : (0.1,  50.0),
    "flow_organic_mL_min"   : (0.1,  20.0),
    "FRR"                   : (1.0,   5.0),
    "mixer_flow_mL_min"     : (0.2,  70.0),
    "dilution_flow_mL_min"  : (0.1, 250.0),
    "tff_feed_mL_min"       : (0.3, 320.0),
    "retentate_mL_min"      : (0.05, 55.0),
    "permeate_mL_min"       : (0.2, 300.0),
    "after_cryo_flow_mL_min": (0.15, 57.0),
    "final_product_mL_min"  : (0.1,  56.0),
    "dilution_ratio"        : (0.5,  10.0),
    "tff_vrr"               : (2.0,  20.0),
    "cryo_flow_mL_min"      : (0.05,  5.0),
    "sterile_filter_dead_vol_mL_min": (0.01, 5.0),
}

# Typical (expected) operating ranges — tighter than the hard limits above
TYPICAL_RANGES = {
    "flow_aqueous_mL_min"   : (3.0,  15.0),
    "flow_organic_mL_min"   : (1.0,   5.0),
    "FRR"                   : (1.0,   5.0),
    "mixer_flow_mL_min"     : (4.0,  20.0),
    "dilution_flow_mL_min"  : (4.0,  75.0),
    "tff_feed_mL_min"       : (8.0, 100.0),
    "retentate_mL_min"      : (1.0,  15.0),
    "permeate_mL_min"       : (5.0,  90.0),
    "after_cryo_flow_mL_min": (1.0,  17.0),
    "final_product_mL_min"  : (1.0,  16.0),
    "dilution_ratio"        : (1.0,   5.0),
    "tff_vrr"               : (4.0,  10.0),
    "cryo_flow_mL_min"      : (0.1,   2.0),
    "sterile_filter_dead_vol_mL_min": (0.1, 1.0),
}

N_FLOW   = 9
N_EQ     = 7
BASE_DOF = N_FLOW - N_EQ   # = 2
TOL      = 1e-6


# ── A: Input Reader ─────────────────────────────────────────────────────────────
# Reads the new compact input format.
# Lines with | separators: only the text before the first | is parsed.
# Blank after = means unknown. Any number means known.

def read_input(filepath):
    """
    Reads the unified input file (lnp_input.txt).
    Parses both Section A (flow rates) and Section B (composition).

    Returns:
      known      : dict of known flow-rate variables  {var_name: float}
      params     : dict of required process parameters {param_name: float}
      comp_known : dict of composition inputs  {(stream, comp, field): float}
      raw        : list of all file lines
    """
    known, params, comp_known, raw = {}, {}, {}, []
    with open(filepath, encoding="utf-8") as f:
        raw = f.readlines()

    valid_comp_streams = {"aqueous", "organic", "dilution_buffer", "cryo_feed", "final_product"}
    valid_comp_fields  = {"conc_mg_mL", "mass_mg_min"}
    valid_comp_names   = {
        "mRNA", "Citric_Acid", "Sodium_Citrate", "EDTA", "Water",
        "Ionisable_Lipid", "DSPC", "Cholesterol", "PEG_Lipid", "Ethanol",
        "NaCl", "Na2HPO4", "NaH2PO4", "Trehalose",
    }

    for line in raw:
        parse_part = line.split("|")[0].strip()

        if not parse_part or "=" not in parse_part:
            continue

        key, _, val = parse_part.partition("=")
        key, val = key.strip(), val.strip()

        if not key:
            continue

        if key in ALL_VARIABLES:
            if val:
                try:
                    known[key] = float(val)
                except ValueError:
                    print(f"  [WARNING] Cannot parse '{key}' = '{val}'. Treating as UNKNOWN.")

        elif key in PROCESS_PARAMS:
            if not val:
                raise ValueError(
                    f"\n  [ERROR] Required parameter '{key}' has no value.\n"
                    f"  This is a fixed operating condition and must always be provided."
                )
            try:
                params[key] = float(val)
            except ValueError:
                raise ValueError(f"\n  [ERROR] Invalid value for '{key}': '{val}'")

        elif "." in key:
            parts = key.split(".")
            if len(parts) == 3:
                stream, comp, field = parts
                if (stream in valid_comp_streams
                        and comp in valid_comp_names
                        and field in valid_comp_fields
                        and val):
                    try:
                        comp_known[(stream, comp, field)] = float(val)
                    except ValueError:
                        pass

    for p in PROCESS_PARAMS:
        if p not in params:
            raise ValueError(
                f"\n  [ERROR] Required parameter '{p}' is missing from the input file."
            )
    return known, params, comp_known, raw


def derive_flows_from_composition(comp_known, params):
    """
    Reverse path: derives stream flow rates from composition data.

    For each primary stream, if the user gave both conc (mg/mL) and mass (mg/min)
    for any component, the implied flow = mass / conc.

    Multiple components in the same stream should all imply the same flow rate.
    The MEDIAN of all implied values is used — this is robust to a single
    outlier (e.g., a typo in one mass value). Any component whose implied flow
    deviates more than 1% from the median is reported as suspect.

    After the flow rate is determined, all mass values are corrected to
    conc x flow in the composition solver. So a wrong mass entry does not
    propagate — it gets overwritten.

    Only runs for streams where the flow rate is not already given in Section A.
    """
    stream_to_flow_var = {
        "aqueous"         : "flow_aqueous_mL_min",
        "organic"         : "flow_organic_mL_min",
        "dilution_buffer" : "dilution_flow_mL_min",
        "cryo_feed"       : "cryo_flow_mL_min",
    }
    stream_comps = {
        "aqueous"         : ["mRNA", "Citric_Acid", "Sodium_Citrate", "EDTA", "Water"],
        "organic"         : ["Ionisable_Lipid", "DSPC", "Cholesterol", "PEG_Lipid", "Ethanol"],
        "dilution_buffer" : ["NaCl", "Na2HPO4", "NaH2PO4", "Water"],
        "cryo_feed"       : ["Trehalose", "Water"],
    }

    derived = {}

    for stream, flow_var in stream_to_flow_var.items():
        implied = []
        for comp in stream_comps[stream]:
            conc = comp_known.get((stream, comp, "conc_mg_mL"))
            mass = comp_known.get((stream, comp, "mass_mg_min"))
            if conc is not None and mass is not None and conc > 0:
                implied.append((comp, mass / conc))

        if not implied:
            continue

        flows_only = sorted(f for _, f in implied)
        median_flow = flows_only[len(flows_only) // 2]

        suspect = [(c, f) for c, f in implied
                   if median_flow > 0 and abs(f - median_flow) / median_flow > 0.01]
        for comp, bad_flow in suspect:
            print(f"  [INFO] Reverse flow derivation: '{stream}.{comp}' implies "
                  f"flow = {bad_flow:.4f} mL/min, which differs from the median "
                  f"({median_flow:.4f} mL/min) by more than 1%. "
                  f"This is likely a wrong mass value. It will be corrected to "
                  f"conc x flow after solving.")

        derived[flow_var] = round(median_flow, 6)

    return derived


# ── B: Physical Validation (Input Values) ──────────────────────────────────────
# Checks user-provided inputs BEFORE solving.
# Returns list of warning dicts.

def validate_inputs(known_vars, params):
    warnings = []

    all_to_check = {**known_vars, **params}

    for name, val in all_to_check.items():
        # Negative value check (physically impossible)
        if val < 0:
            warnings.append({
                "level"  : "ERROR",
                "variable": name,
                "value"  : val,
                "message": f"NEGATIVE VALUE detected for '{DISPLAY.get(name, name)}' "
                           f"(value = {val}). Flow rates cannot be negative. "
                           f"Please correct this in the input file."
            })
            continue

        # Hard realistic range check
        if name in REALISTIC_RANGES:
            rmin, rmax = REALISTIC_RANGES[name]
            if val < rmin or val > rmax:
                warnings.append({
                    "level"  : "WARNING",
                    "variable": name,
                    "value"  : val,
                    "message": f"Value {val} for '{DISPLAY.get(name, name)}' is outside "
                               f"the realistic operating range [{rmin}, {rmax}]. "
                               f"This may cause physically meaningless results."
                })
                continue  # Don't also flag as typical-range if already outside realistic

        # Typical range advisory
        if name in TYPICAL_RANGES:
            tmin, tmax = TYPICAL_RANGES[name]
            if not (tmin <= val <= tmax):
                rmin, rmax = REALISTIC_RANGES.get(name, (None, None))
                if rmin is not None and rmin <= val <= rmax:
                    warnings.append({
                        "level"  : "ADVISORY",
                        "variable": name,
                        "value"  : val,
                        "message": f"Value {val} for '{DISPLAY.get(name, name)}' is "
                                   f"outside the typical operating range [{tmin}, {tmax}] "
                                   f"but within realistic limits. Proceed with caution."
                    })

    return warnings


# ── C: Physical Validation (Solved Results) ────────────────────────────────────
# Checks ALL solved variables after the calculation is complete.

def validate_solution(solution, params):
    warnings = []

    # Check every solved flow variable
    for name in ALL_VARIABLES:
        val = solution.get(name)
        if val is None:
            continue

        # Negative flow: physically impossible
        if val < 0:
            warnings.append({
                "level"  : "ERROR",
                "variable": name,
                "value"  : val,
                "message": f"NEGATIVE FLOW: '{DISPLAY.get(name, name)}' = {val:.4f} mL/min. "
                           f"This is physically impossible. Your input values are likely "
                           f"inconsistent or outside realistic ranges."
            })
            continue

        # Hard realistic range
        if name in REALISTIC_RANGES:
            rmin, rmax = REALISTIC_RANGES[name]
            if val < rmin or val > rmax:
                warnings.append({
                    "level"  : "WARNING",
                    "variable": name,
                    "value"  : val,
                    "message": f"'{DISPLAY.get(name, name)}' = {val:.4f} mL/min is outside "
                               f"the realistic operating range [{rmin}, {rmax}]. "
                               f"Review your inputs."
                })
                continue

        # Typical range advisory
        if name in TYPICAL_RANGES:
            tmin, tmax = TYPICAL_RANGES[name]
            if not (tmin <= val <= tmax):
                warnings.append({
                    "level"  : "ADVISORY",
                    "variable": name,
                    "value"  : val,
                    "message": f"'{DISPLAY.get(name, name)}' = {val:.4f} mL/min is outside "
                               f"the typical range [{tmin}, {tmax}]. Acceptable but unusual."
                })

    # FRR specific check
    frr = solution.get("FRR", float("nan"))
    if not (1.0 <= frr <= 5.0):
        warnings.append({
            "level"  : "WARNING",
            "variable": "FRR",
            "value"  : frr,
            "message": f"FRR = {frr:.4f} is outside the typical range [1, 5]. "
                       f"Very low or very high FRR can lead to poor LNP formation "
                       f"or unstable particle size distribution."
        })

    # Retentate must be less than TFF feed
    ret = solution.get("retentate_mL_min", 0)
    tff = solution.get("tff_feed_mL_min", 1)
    if ret >= tff:
        warnings.append({
            "level"  : "ERROR",
            "variable": "retentate_mL_min",
            "value"  : ret,
            "message": f"Retentate ({ret:.4f}) >= TFF Feed ({tff:.4f}). "
                       f"Retentate must be smaller than TFF Feed. "
                       f"Check tff_vrr: it must be > 1."
        })

    # Permeate must be positive
    perm = solution.get("permeate_mL_min", 0)
    if perm <= 0:
        warnings.append({
            "level"  : "ERROR",
            "variable": "permeate_mL_min",
            "value"  : perm,
            "message": f"Permeate = {perm:.4f}. TFF permeate must be positive. "
                       f"Check tff_vrr value — VRR must be > 1."
        })

    # Low process yield advisory
    aq   = solution.get("flow_aqueous_mL_min", 0)
    org  = solution.get("flow_organic_mL_min", 0)
    dil  = solution.get("dilution_flow_mL_min", 0)
    cryo = params.get("cryo_flow_mL_min", 0)
    prod = solution.get("final_product_mL_min", 0)
    total_in = aq + org + dil + cryo
    if total_in > 0:
        yield_pct = 100.0 * prod / total_in
        if yield_pct < 5.0:
            warnings.append({
                "level"  : "ADVISORY",
                "variable": "final_product_mL_min",
                "value"  : prod,
                "message": f"Final product flow is very low relative to total process input. "
                           f"This may be expected for high-VRR runs, but verify your "
                           f"tff_vrr and sterile_filter_dead_vol values."
            })

    return warnings


# ── D: FRR Resolution ───────────────────────────────────────────────────────────

def resolve_frr(kv, log):
    aq  = kv.get("flow_aqueous_mL_min")
    org = kv.get("flow_organic_mL_min")
    frr = kv.get("FRR")
    mix = kv.get("mixer_flow_mL_min")

    if aq is not None and org is not None and org != 0:
        c = aq / org
        if frr is None:
            kv["FRR"] = c
            log.append(f"    FRR derived: Aq/Org = {aq}/{org} = {c:.4f}")
        else:
            log.append(f"    FRR check: given={frr}, Aq/Org={c:.4f}, diff={abs(c-frr):.2e}")
    elif frr is not None and org is not None and aq is None:
        kv["flow_aqueous_mL_min"] = frr * org
        log.append(f"    Aq derived: FRR x Org = {frr}x{org} = {frr*org:.4f}")
    elif frr is not None and aq is not None and org is None:
        kv["flow_organic_mL_min"] = aq / frr
        log.append(f"    Org derived: Aq/FRR = {aq}/{frr} = {aq/frr:.4f}")
    elif frr is not None and mix is not None and aq is None and org is None:
        kv["flow_aqueous_mL_min"] = frr / (1 + frr) * mix
        kv["flow_organic_mL_min"] = mix / (1 + frr)
        log.append(f"    Aq,Org from FRR+Mix: Aq={kv['flow_aqueous_mL_min']:.4f}, "
                   f"Org={kv['flow_organic_mL_min']:.4f}")
    else:
        log.append("    FRR: will be computed post-solve from Aq/Org.")


# ── E: DOF Analysis ─────────────────────────────────────────────────────────────

def dof_analysis(kv):
    n_known = sum(1 for v in FLOW_VARIABLES if v in kv)
    return BASE_DOF - n_known, n_known


# ── F: Redundant Variable Selector ──────────────────────────────────────────────

def pick_redundant(kv, n, log):
    cands = sorted(
        [(k, KNOWABILITY.get(k, 0.5)) for k in kv if k in FLOW_VARIABLES],
        key=lambda x: x[1]
    )
    chosen = [name for name, _ in cands[:n]]
    for name, sc in cands[:n]:
        log.append(f"    Assumed unknown: '{DISPLAY.get(name, name)}' (knowability = {sc:.2f})")
    return chosen


# ── G: Equation System Builder ──────────────────────────────────────────────────
# 7 equations x 9 variables
# Vars: [Aq, Org, Mix, Dil, TFF_in, Ret, Perm, Cryo_out, Prod]
#         0    1    2    3     4      5    6       7         8

def build_system(params):
    d    = params["dilution_ratio"]
    vrr  = params["tff_vrr"]
    cryo = params["cryo_flow_mL_min"]
    dead = params["sterile_filter_dead_vol_mL_min"]

    A = np.zeros((7, 9));  b = np.zeros(7)

    A[0, 0]=1;  A[0, 1]=1;  A[0, 2]=-1          # EQ1: Aq+Org=Mix
    A[1, 2]=d;  A[1, 3]=-1                        # EQ2: d*Mix=Dil
    A[2, 2]=1;  A[2, 3]=1;  A[2, 4]=-1           # EQ3: Mix+Dil=TFF
    A[3, 4]=1;  A[3, 5]=-vrr                      # EQ4: TFF=vrr*Ret
    A[4, 4]=1;  A[4, 5]=-1; A[4, 6]=-1           # EQ5: TFF=Ret+Perm
    A[5, 5]=1;  A[5, 7]=-1; b[5]=-cryo           # EQ6: Ret+cryo=Cryo_out
    A[6, 7]=1;  A[6, 8]=-1; b[6]=dead            # EQ7: Cryo_out-dead=Prod

    return A, b


# ── H: Unified Solver ───────────────────────────────────────────────────────────

def solve(known_vars_input, params):
    log = []
    assumed_unknown = []
    verification = {}

    kv = dict(known_vars_input)

    log.append("  [Step 1] FRR Resolution:")
    resolve_frr(kv, log)

    known_flow = {k: v for k, v in kv.items() if k in FLOW_VARIABLES}
    dof, n_known = dof_analysis(known_flow)

    log.append(f"\n  [Step 2] DOF Analysis:")
    log.append(f"    Flow variables : {N_FLOW}  |  Equations : {N_EQ}  |  Base DOF : {BASE_DOF}")
    log.append(f"    Known flow vars: {n_known}  |  Effective DOF : {dof}")

    if dof > 0:
        unknown_list = [v for v in FLOW_VARIABLES if v not in known_flow]
        raise ValueError(
            f"\n  [SOLVER ERROR] UNDER-DETERMINED SYSTEM  (DOF = {dof})\n"
            f"  Need {dof} more known flow value(s).\n"
            f"  Currently unknown: {unknown_list}\n"
            f"  Please fill in {dof} more value(s) in lnp_input.txt and re-run."
        )

    if dof < 0:
        n_red = abs(dof)
        log.append(f"\n  [Step 3] Over-determined by {n_red}. Selecting redundant variable(s)...")
        assumed_unknown = pick_redundant(known_flow, n_red, log)
        verify_targets  = {k: kv[k] for k in assumed_unknown}
        for k in assumed_unknown:
            del kv[k]
        known_flow = {k: v for k, v in kv.items() if k in FLOW_VARIABLES}
    else:
        log.append(f"\n  [Step 3] Exactly determined. No redundancy handling needed.")
        verify_targets = {}

    A_full, b_full = build_system(params)
    var_list  = FLOW_VARIABLES
    known_idx = {var_list.index(k): v for k, v in known_flow.items()}
    unk_idx   = [i for i in range(N_FLOW) if i not in known_idx]
    unk_names = [var_list[i] for i in unk_idx]

    b_mod = b_full.copy()
    for idx, val in known_idx.items():
        b_mod -= A_full[:, idx] * val
    A_red = A_full[:, unk_idx]

    log.append(f"\n  [Step 4] Linear system solve:")
    log.append(f"    Known   : {list(known_flow.keys())}")
    log.append(f"    Unknown : {unk_names}")
    log.append(f"    Solver  : numpy.linalg.lstsq (linear least-squares)")
    log.append(f"    All equations are LINEAR — no iterative method needed.")

    if len(unk_idx) == 0:
        x_sol = np.array([])
    else:
        x_sol, residuals, rank, sv = np.linalg.lstsq(A_red, b_mod, rcond=None)
        cond = sv[0]/sv[-1] if len(sv)>1 and sv[-1]!=0 else float("inf")
        res_ok = (len(residuals)==0 or np.all(np.array(residuals)<1e-8))
        log.append(f"    Rank: {rank}  |  Condition number: {cond:.2f}  |  Residuals: {'Negligible' if res_ok else residuals}")

    # Assemble solution
    solution = {k: v for k, v in kv.items() if k in FLOW_VARIABLES}
    for i, idx in enumerate(unk_idx):
        solution[var_list[idx]] = float(x_sol[i])

    # Compute FRR post-solve
    aq  = solution.get("flow_aqueous_mL_min", 0.0)
    org = solution.get("flow_organic_mL_min", 0.0)
    solution["FRR"] = aq / org if org != 0 else kv.get("FRR", float("nan"))
    log.append(f"\n  [Step 5] FRR post-solve: {aq:.4f} / {org:.4f} = {solution['FRR']:.4f}")

    # Redundancy verification
    if assumed_unknown:
        log.append(f"\n  [Step 6] Redundancy verification:")
        all_pass = True
        for vname in assumed_unknown:
            user_val   = verify_targets[vname]
            solved_val = solution.get(vname, float("nan"))
            diff       = abs(solved_val - user_val)
            passed     = diff < TOL
            if not passed: all_pass = False
            status = "PASS" if passed else f"FAIL (diff={diff:.2e})"
            verification[vname] = {
                "user_provided": user_val, "solver_result": solved_val,
                "difference": diff, "passed": passed, "status": status
            }
            log.append(f"    {DISPLAY.get(vname, vname)}: user={user_val:.6f}, "
                       f"solver={solved_val:.6f} -> {status}")
        log.append(f"    Result: {'Consistent' if all_pass else 'INCONSISTENCY DETECTED — review inputs.'}")

    return solution, assumed_unknown, verification, log, dof


# ── I: Material Balance Verifier ────────────────────────────────────────────────

def verify_balance(sol, params):
    d    = params["dilution_ratio"]
    vrr  = params["tff_vrr"]
    cryo = params["cryo_flow_mL_min"]
    dead = params["sterile_filter_dead_vol_mL_min"]

    aq=sol["flow_aqueous_mL_min"]; org=sol["flow_organic_mL_min"]
    frr=sol["FRR"]; mix=sol["mixer_flow_mL_min"]; dil=sol["dilution_flow_mL_min"]
    tff=sol["tff_feed_mL_min"];  ret=sol["retentate_mL_min"]; perm=sol["permeate_mL_min"]
    cout=sol["after_cryo_flow_mL_min"]; prod=sol["final_product_mL_min"]

    checks = []
    def chk(name, desc, lhs, rhs):
        err = abs(lhs - rhs)
        checks.append({"name":name,"description":desc,"LHS":lhs,"RHS":rhs,"error":err,"passed":err<TOL})

    chk("EQ1",    "Aq + Org = Mixer Output",              aq+org,      mix)
    chk("EQ-FRR", "Aq = FRR x Org",                      aq,          frr*org)
    chk("EQ2",    f"Dil_ratio({d}) x Mix = Dil",         d*mix,       dil)
    chk("EQ3",    "Mix + Dil = TFF Feed",                 mix+dil,     tff)
    chk("EQ4",    f"TFF Feed = VRR({vrr}) x Retentate",  tff,         vrr*ret)
    chk("EQ5",    "TFF Feed = Retentate + Permeate",      tff,         ret+perm)
    chk("EQ6",    f"Ret + {cryo}(cryo) = After-Cryo",    ret+cryo,    cout)
    chk("EQ7",    f"After-Cryo - {dead}(dead) = Product", cout-dead,  prod)
    chk("OVERALL","Total In = Total Out",                 aq+org+dil+cryo, prod+perm+dead)
    return checks


# ── J: Calculation Mode Identifier ──────────────────────────────────────────────

def identify_mode(kvi, au):
    eff = set(kvi) - set(au)
    if {"flow_aqueous_mL_min","flow_organic_mL_min"}.issubset(eff) and "final_product_mL_min" not in eff:
        return "FORWARD", ("FORWARD CALCULATION\n"
            "  Upstream pump flows specified. Simulation propagates forward\n"
            "  through the process to calculate all downstream variables.")
    if "final_product_mL_min" in eff and not {"flow_aqueous_mL_min","flow_organic_mL_min"}.issubset(eff):
        return "BACKWARD", ("BACKWARD CALCULATION\n"
            "  Final Product flow specified. Simulation works backwards to\n"
            "  determine the required upstream flow rates.")
    return "GENERAL", ("GENERAL / MID-PROCESS CALCULATION\n"
        "  A mix of variables from different process sections specified.\n"
        "  Simulation solves material balance equations simultaneously.")


# ── K: Output Folder ─────────────────────────────────────────────────────────────

def create_folder():
    base = "simulation_runs"
    os.makedirs(base, exist_ok=True)
    run_id = len(os.listdir(base)) + 1
    folder = os.path.join(base, f"run_{run_id:03d}")
    os.makedirs(folder, exist_ok=True)
    return folder, run_id


# ── L: Warning Formatter ─────────────────────────────────────────────────────────

LEVEL_SYMBOL = {"ERROR": "!! ERROR  !!", "WARNING": "** WARNING **", "ADVISORY": "-- ADVISORY --"}

def format_warnings_block(warnings, title):
    if not warnings:
        return []
    L = []
    L.append(f"  {title}")
    errors     = [w for w in warnings if w["level"] == "ERROR"]
    warns      = [w for w in warnings if w["level"] == "WARNING"]
    advisories = [w for w in warnings if w["level"] == "ADVISORY"]
    for group in [errors, warns, advisories]:
        for w in group:
            sym = LEVEL_SYMBOL.get(w["level"], w["level"])
            L.append(f"  [{sym}]  {w['message']}")
    L.append("")
    return L


# ── M: Output Report Writer ──────────────────────────────────────────────────────

def write_report(folder, run_id, sol, params, kvi, au, verif, dof,
                 mode, mdesc, bchk, input_warnings, solution_warnings):
    now = datetime.now().strftime("%Y-%m-%d  %H:%M:%S")
    L=[]
    def a(s=""): L.append(s)

    a("  LNP MANUFACTURING PROCESS SIMULATION — OUTPUT REPORT")
    a(f"  Run ID : {run_id:03d}    |    Timestamp : {now}")
    a()

    # Calculation mode
    a("  CALCULATION MODE:")
    for ln in mdesc.split("\n"): a("  "+ln)
    a()

    # Input warnings
    if input_warnings:
        for ln in format_warnings_block(input_warnings, "INPUT VALIDATION WARNINGS"):
            a(ln)

    # Solution warnings
    if solution_warnings:
        for ln in format_warnings_block(solution_warnings, "SOLUTION VALIDATION WARNINGS"):
            a(ln)

    # Fixed parameters
    a("  FIXED PROCESS PARAMETERS")
    a(f"  {'Parameter':<44} {'Value':>10}  Unit")
    for k, (disp, unit) in {
        "dilution_ratio":("Dilution Ratio","—"),
        "tff_vrr":("TFF Volume Reduction Ratio (VRR)","—"),
        "cryo_flow_mL_min":("Cryoprotectant Addition Flow Rate","mL/min"),
        "sterile_filter_dead_vol_mL_min":("Sterile Filter Dead Volume","mL/min"),
    }.items():
        a(f"  {disp:<44} {params[k]:>10.4f}  {unit}")
    a()

    # DOF analysis
    a("  DEGREES OF FREEDOM (DOF) ANALYSIS")
    a(f"  Total flow variables        : {N_FLOW}")
    a(f"  Independent equations       : {N_EQ}")
    a(f"  Base DOF                    : {BASE_DOF}")
    a(f"  Known values (user)         : {sum(1 for v in FLOW_VARIABLES if v in kvi)}")
    a(f"  FRR specified by user       : {'Yes' if 'FRR' in kvi else 'No'}")
    a(f"  Assumed unknown (redundancy): {len(au)}")
    a(f"  Effective DOF               : {dof}")
    if dof==0: a("  Status: EXACTLY DETERMINED — unique solution obtained.")
    elif dof<0: a(f"  Status: OVER-DETERMINED by {abs(dof)} — redundancy check performed.")
    a()

    # Input summary
    a("  INPUT SUMMARY")
    a(f"  {'Variable':<47} {'Status':<14} {'User Value':>12}")
    for vn in ALL_VARIABLES:
        disp=DISPLAY.get(vn,vn); val=kvi.get(vn)
        if vn in au:          st="ASSUME-UNKN"; vs=f"{val:.4f}" if val else "—"
        elif val is not None: st="KNOWN";        vs=f"{val:.4f}"
        else:                 st="UNKNOWN";      vs="(to solve)"
        a(f"  {disp:<47} {st:<14} {vs:>12}")
    a()

    # Results by section
    a("  SIMULATION RESULTS — BY PROCESS SECTION")
    a()
    UNITS = {k:"mL/min" for k in FLOW_VARIABLES}; UNITS["FRR"]="dimensionless"
    secs = [
        ("SECTION 1 — MICROFLUIDIC MIXING  (M-101 T-Mixer)",
         ["flow_aqueous_mL_min","flow_organic_mL_min","FRR","mixer_flow_mL_min"]),
        ("SECTION 2 — BUFFER DILUTION  (T-103 Post-Mix Dilution Vessel)",
         ["dilution_flow_mL_min","tff_feed_mL_min"]),
        ("SECTION 3 — TANGENTIAL FLOW FILTRATION  (UF-101, 100 kDa MWCO)",
         ["retentate_mL_min","permeate_mL_min"]),
        ("SECTION 4 — CRYOPROTECTANT ADDITION  (T-100 Stabilizer Vessel)",
         ["after_cryo_flow_mL_min"]),
        ("SECTION 5 — STERILE FILTRATION  (F-101, 0.22 um Filter)",
         ["final_product_mL_min"]),
    ]
    for title, vnames in secs:
        a(f"  {title}")
        a(f"  {'Variable':<47} {'Value':>12}  {'Unit':<16} Source")
        for vn in vnames:
            disp=DISPLAY.get(vn,vn); val=sol.get(vn,float("nan"))
            unit=UNITS.get(vn,"")
            src="(user input)" if (vn in kvi and vn not in au) else "(calculated)"
            a(f"  {disp:<47} {val:>12.4f}  {unit:<16} {src}")
        a()

    # Overall balance
    a("  OVERALL PROCESS MATERIAL BALANCE")
    aq=sol["flow_aqueous_mL_min"]; org=sol["flow_organic_mL_min"]
    dil=sol["dilution_flow_mL_min"]; cryo=params["cryo_flow_mL_min"]
    prod=sol["final_product_mL_min"]; perm=sol["permeate_mL_min"]
    dead=params["sterile_filter_dead_vol_mL_min"]
    ti=aq+org+dil+cryo; to=prod+perm+dead
    a("  INPUTS:")
    a(f"    Aqueous Phase Flow           : {aq:>10.4f}  mL/min")
    a(f"    Organic Phase Flow           : {org:>10.4f}  mL/min")
    a(f"    Dilution Buffer Flow         : {dil:>10.4f}  mL/min")
    a(f"    Cryoprotectant Flow          : {cryo:>10.4f}  mL/min")
    a(f"    TOTAL INPUTS                 : {ti:>10.4f}  mL/min")
    a()
    a("  OUTPUTS:")
    a(f"    Final Product                : {prod:>10.4f}  mL/min")
    a(f"    TFF Permeate (waste)         : {perm:>10.4f}  mL/min")
    a(f"    Sterile Filter Dead Volume   : {dead:>10.4f}  mL/min")
    a(f"    TOTAL OUTPUTS                : {to:>10.4f}  mL/min")
    a()
    imb=abs(ti-to)
    a(f"  Balance Status : {'BALANCED' if imb<TOL else f'IMBALANCED (error={imb:.2e})'}")
    a()

    # Redundancy verification
    if verif:
        a("  REDUNDANCY SELF-VERIFICATION REPORT")
        a("  Variable(s) treated as unknown during solve, then compared to user values:")
        a()
        all_p = all(v["passed"] for v in verif.values())
        for vn, vd in verif.items():
            a(f"  Variable : {DISPLAY.get(vn, vn)}")
            a(f"    User provided value  : {vd['user_provided']:.6f}  mL/min")
            a(f"    Solver calculated    : {vd['solver_result']:.6f}  mL/min")
            a(f"    Difference           : {vd['difference']:.2e}")
            a(f"    Result               : {vd['status']}"); a()
        a(f"  CONCLUSION: {'Data is INTERNALLY CONSISTENT' if all_p else 'INCONSISTENCY FOUND — review inputs.'}")
        a()

    a("  END OF REPORT")

    with open(os.path.join(folder,"output_report.txt"),"w",encoding="utf-8") as f:
        f.write("\n".join(L))


# ── N: Diagnostics Writer ────────────────────────────────────────────────────────

def write_diagnostics(folder, solve_log, bchk, params, dof, mode, au, verif,
                      input_warnings, solution_warnings):
    L=[]
    def a(s=""): L.append(s)

    a("  LNP PROCESS SIMULATION — DIAGNOSTICS REPORT")
    a()
    a("  SOLVER INFORMATION:")
    a("  Library   : numpy.linalg.lstsq (linear least-squares)")
    a("  Method    : Direct matrix decomposition")
    a("  Eq. type  : All equations are LINEAR material balances")
    a("  FRR       : Resolved analytically before matrix build (not in matrix)")
    a("  Iterative : NOT REQUIRED")
    a()

    d=params["dilution_ratio"]; vrr=params["tff_vrr"]
    cryo=params["cryo_flow_mL_min"]; dead=params["sterile_filter_dead_vol_mL_min"]

    a("  GOVERNING EQUATIONS")
    a(f"  {'ID':<10} {'Name':<30} Equation")
    for row in [
        ("EQ1",    "Mixer Balance",        "Aq + Org = Mix"),
        ("EQ-FRR", "FRR (analytical)",     "FRR = Aq / Org"),
        ("EQ2",    "Dilution Flow",        f"Dil = {d} x Mix"),
        ("EQ3",    "Dilution Balance",     "Mix + Dil = TFF_Feed"),
        ("EQ4",    "TFF VRR",             f"TFF_Feed = {vrr} x Retentate"),
        ("EQ5",    "TFF Balance",          "TFF_Feed = Retentate + Permeate"),
        ("EQ6",    "Cryo Addition",       f"Ret + {cryo} = After_Cryo"),
        ("EQ7",    "Sterile Filter",      f"After_Cryo - {dead} = Final_Product"),
        ("OVERALL","Overall Balance",      "Aq+Org+Dil+Cryo = Prod+Perm+Dead"),
    ]: a(f"  {row[0]:<10} {row[1]:<30} {row[2]}")
    a()

    a("  EQUATIONS PER PROCESS SECTION")
    for sec, desc in [
        ("Mixing (M-101)",         "2 eqs -> EQ1: Aq+Org=Mix  |  EQ-FRR: FRR=Aq/Org"),
        ("Dilution (T-103)",       "2 eqs -> EQ2: Dil=d*Mix   |  EQ3: Mix+Dil=TFF_Feed"),
        ("TFF (UF-101)",           "2 eqs -> EQ4: TFF=vrr*Ret |  EQ5: TFF=Ret+Perm"),
        ("Cryo Addition (T-100)",  "1 eq  -> EQ6: Ret+cryo=After_Cryo"),
        ("Sterile Filter (F-101)", "1 eq  -> EQ7: After_Cryo-dead=Product"),
        ("Overall Process",        "1 eq  -> Total In = Total Out"),
    ]: a(f"  {sec:<28} : {desc}")
    a()

    # Warnings section
    all_warnings = input_warnings + solution_warnings
    if all_warnings:
        a("  ALL VALIDATION WARNINGS")
        for w in all_warnings:
            sym = LEVEL_SYMBOL.get(w["level"], w["level"])
            a(f"  [{sym}]  {w['message']}")
        a()

    a("  SOLVE LOG")
    for ln in solve_log: a(ln)
    a()

    a("  MATERIAL BALANCE VERIFICATION")
    a(f"  {'Eq':<10} {'Description':<42} {'LHS':>10} {'RHS':>10} {'Error':>12}  Status")
    for c in bchk:
        sym="PASS" if c["passed"] else "FAIL"
        a(f"  {c['name']:<10} {c['description']:<42} {c['LHS']:>10.4f} "
          f"{c['RHS']:>10.4f} {c['error']:>12.2e}  {sym}")
    a()
    a(f"  RESULT: {'ALL SATISFIED' if all(c['passed'] for c in bchk) else 'SOME FAILED'}")
    a()

    if au:
        a("  KNOWABILITY ANALYSIS — REDUNDANCY RESOLUTION")
        a("  Score: 0.0 = least independently known | 1.0 = set directly by user")
        a()
        for vn in au:
            sc=KNOWABILITY.get(vn,0.5); disp=DISPLAY.get(vn,vn); vd=verif.get(vn,{})
            if sc<0.55:   reason="Typically derived; rarely directly measured."
            elif sc<0.75: reason="Can be measured but often computed from other variables."
            else:         reason="Generally well-known; chosen to balance DOF."
            a(f"  Variable : {disp}")
            a(f"    Knowability : {sc:.2f}  |  Reason: {reason}")
            if vd:
                a(f"    User value  : {vd['user_provided']:.6f}  |  Solver: {vd['solver_result']:.6f}  |  {vd['status']}")
            a()

    a("  END OF DIAGNOSTICS")
    with open(os.path.join(folder,"diagnostics.txt"),"w",encoding="utf-8") as f:
        f.write("\n".join(L))


# ── O: Main Entry Point ──────────────────────────────────────────────────────────

def run_simulation(input_file=INPUT_FILE):
    print()
    print("  LNP MANUFACTURING PROCESS SIMULATION")
    print(f"  Input file   : {input_file}")
    print(f"  User manual  : lnp_user_manual.txt")
    print()

    try:
        kvi, params, comp_known, raw = read_input(input_file)
    except FileNotFoundError:
        print(f"  [ERROR] '{input_file}' not found.")
        print("  Place lnp_input.txt in the same folder as this script and re-run.")
        return None
    except ValueError as e:
        print(str(e)); return None

    if comp_known:
        derived_flows = derive_flows_from_composition(comp_known, params)
        for fv, fval in derived_flows.items():
            if fv not in kvi:
                kvi[fv] = fval
                print(f"  [INFO] Flow rate derived from composition: {DISPLAY.get(fv, fv)} = {fval:.4f} mL/min")
        if derived_flows:
            print()

    input_warnings = validate_inputs(kvi, params)
    errors_in_input = [w for w in input_warnings if w["level"] == "ERROR"]

    print(f"  Variables provided : {len(kvi)}")
    for k, v in kvi.items():
        print(f"    {DISPLAY.get(k,k)} = {v}")
    print()

    if input_warnings:
        print(f"  INPUT VALIDATION ({len(input_warnings)} message(s)):")
        for w in input_warnings:
            sym = LEVEL_SYMBOL.get(w["level"], w["level"])
            print(f"    [{sym}] {w['message']}")
        print()

    if errors_in_input:
        print("  [STOPPED] Critical input errors found. Please fix them before re-running.")
        return None

    try:
        sol, au, verif, solve_log, dof = solve(kvi, params)
    except ValueError as e:
        print(str(e)); return None

    solution_warnings = validate_solution(sol, params)

    mode, mdesc = identify_mode(kvi, au)
    bchk = verify_balance(sol, params)
    all_bal = all(c["passed"] for c in bchk)

    folder, run_id = create_folder()
    shutil.copy(input_file, os.path.join(folder,"input_used.txt"))
    write_report(folder, run_id, sol, params, kvi, au, verif, dof,
                 mode, mdesc, bchk, input_warnings, solution_warnings)
    write_diagnostics(folder, solve_log, bchk, params, dof, mode, au, verif,
                      input_warnings, solution_warnings)

    print(f"  Calculation Mode : {mode}  |  DOF : {dof}")
    print()
    print("  RESULTS:")
    for vn in ALL_VARIABLES:
        disp=DISPLAY.get(vn,vn); val=sol.get(vn,float("nan"))
        src=" (user input)" if (vn in kvi and vn not in au) else ""
        print(f"    {disp:<47} : {val:>10.4f}{src}")

    aq=sol["flow_aqueous_mL_min"]; org=sol["flow_organic_mL_min"]
    dil=sol["dilution_flow_mL_min"]; cryo=params["cryo_flow_mL_min"]
    prod=sol["final_product_mL_min"]; perm=sol["permeate_mL_min"]
    dead=params["sterile_filter_dead_vol_mL_min"]
    print()
    print(f"  Overall Balance  : {'BALANCED' if all_bal else 'CHECK DIAGNOSTICS'}")

    if solution_warnings:
        print()
        print(f"  SOLUTION WARNINGS ({len(solution_warnings)} message(s)):")
        for w in solution_warnings:
            sym = LEVEL_SYMBOL.get(w["level"], w["level"])
            print(f"    [{sym}] {w['message']}")

    if au:
        all_v = all(v["passed"] for v in verif.values())
        print(f"\n  Redundancy Check : {'ALL PASS' if all_v else 'SOME FAILED — review inputs'}")

    print()
    print(f"  Results saved in : {folder}/")
    print(f"    output_report.txt  |  diagnostics.txt  |  input_used.txt")
    print()

    if comp_known:
        print("  Running composition & mass calculations...")
        print()
        try:
            from lnp_composition import run_composition_from_solver
            run_composition_from_solver(comp_known, sol, params, folder, run_id)
        except ImportError:
            print("  [INFO] lnp_composition.py not found. Skipping composition calculations.")
    else:
        print("  [INFO] No composition data found in Sections B or C of the input file.")
        print("         Fill in concentrations or mass flows there to get composition results.")

    print()
    print("  Simulation complete.")
    print()
    return folder


def zip_all_runs():
    shutil.make_archive("LNP_simulation_results", "zip", "simulation_runs")
    print("  ZIP created: LNP_simulation_results.zip")


if __name__ == "__main__":
    run_simulation(INPUT_FILE)
    # zip_all_runs()  # Uncomment to ZIP all runs
