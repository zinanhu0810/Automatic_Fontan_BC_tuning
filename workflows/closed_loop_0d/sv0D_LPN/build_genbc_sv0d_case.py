#!/usr/bin/env python3
"""Build a svZeroDSolver case mirroring the genBC Fontan closed-loop LPN.

The vascular network is translated directly from the equations and final
parameters in passive_fluid_truncate/genBC_QP_Baretta. The closed-loop SVC/IVC
equations from the legacy tuning workflow are retained in this builder.
Pressures are converted from mmHg to cgs by svZeroDSolver convention.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "sv0D_LPN"
JSON_PATH = OUT_DIR / "Fontan_genBC_1to1_sv0D.json"
NOTEBOOK_PATH = OUT_DIR / "compare_genBC_sv0D.ipynb"

P_CONV = 1334.0


PARAMS = {
    "BSA": 0.85,
    "HR": 90.0,
    "tsvs": 0.1013,
    "tsas": 0.1947,
    "t1": 0.0714,
    "a": 0.137,
    "b": 0.914,
    "csv": 0.102,
    "dsv": 0.074,
    "Vsv0": 18.7802,
    "Ccsa": 6.607,
    "csa": 0.035,
    "dsa": 0.062,
    "Vsa0": 0.5148,
    "P0ith": 0.0,
    "Apith": 0.0,
    "Trr": 0.67,
    "pi": 3.141593,
    "fc": 131.58,
    "fr": 32.895,
    "Rmyo": 0.0085,
    "Lav": 1.66e-05,
    "Kav": 1.88e-04,
    "Kao": 1.97e-05,
    "Cao": 0.6,
    "Ruba": 0.761,
    "Luba": 6.06e-04,
    "Caub": 0.236,
    "Rub": 3.193,
    "Cvub": 3.1472,
    "Rubv": 0.057,
    "Csvc": 1.0277,
    "Rsvc": 0.0076,
    "Rthao": 0.047,
    "Lthao": 0.0018,
    "Cthao": 0.1773,
    "Rabao": 0.646,
    "Labao": 0.0018,
    "Cabao": 0.3,
    "Rlega": 0.929,
    "Llega": 0.0018,
    "Clega": 1.5829,
    "Rlegc": 2.36,
    "Clegv": 7.2032,
    "Rlegv": 0.548,
    "Cabivc": 6.502,
    "Rabivc": 0.0031,
    "V0abivc": 11.1282,
    "R0abivc": 0.04,
    "P0iab": 0.0,
    "Apiab": 0.0,
    "Cthivc": 0.72712,
    "Rthivc": 0.003,
    "V0thivc": 7.2898,
    "R0thivc": 0.01,
    "Llc": 0.00516,
    "Luc": 0.00516,
    "Klc": 0.000516,
    "Kuc": 0.000516,
    "Rla": 13.257,
    "Cl": 4.42,
    "Rlv": 0.045,
    "Rka": 9.201,
    "Ck": 1.83,
    "Rkv": 0.878,
    "Ria": 18.599,
    "Ci": 0.99,
    "Riv": 0.143,
    "Ctcpc": 0.5,
    "Rapoll": 0.073,
    "Rapolr": 0.02057,
    "Lapoll": 0.00227,
    "Lapolr": 0.00227,
    "Clungl": 0.8335,
    "Clungr": 1.0043,
    "Rlungl": 0.121,
    "Rlungr": 0.05569,
    "Llungl": 0.00153,
    "Llungr": 0.00153,
    "Cvpl": 1.5904,
    "Cvpr": 1.8495,
    "Rvpl": 0.034,
    "Rvpr": 0.01439,
    "Clpa": 1.3111,
    "Crpa": 1.3111,
    "coll": 1.0,
    "Rcoll": 5.636,
    "coll_split": 0.5,
    "Cub": 4.184140958,
    "fS": 0.797,
    "fD": 1.0,
    "fCAo": 1.21,
}


X0_MMHG = [
    84.91167,
    39.37003,
    55.49430,
    12.3695,
    12.39335,
    27.36710,
    74.04834,
    24.243500,
    55.40978,
    10.503520,
    45.0045,
    18.56190,
    12.61708,
    12.59937,
    12.06896,
    11.8046,
    16.02366,
    12.25075,
    11.809034,
    11.05807,
    10.000172,
    9.849035,
    10.724361,
    11.49486,
    14.955302,
    10.07034,
    72.11296,
    0.0,
    43.5820,
    -26.0,
    -15.0,
    0.0,
]


def r(name: str) -> float:
    return PARAMS[name] * P_CONV


def c(name: str) -> float:
    return PARAMS[name] / P_CONV


def l(name: str) -> float:
    return PARAMS[name] * P_CONV


def vessel(vid: int, name: str, values: dict[str, float]) -> dict:
    return {
        "vessel_id": vid,
        "vessel_length": 1.0,
        "vessel_name": name,
        "zero_d_element_type": "BloodVessel",
        "zero_d_element_values": values,
    }


def build_config() -> dict:
    p = PARAMS
    split = p["coll_split"]
    rcoll_l = p["Rcoll"] / max(1.0 - split, 1.0e-12)
    rcoll_r = p["Rcoll"] / max(split, 1.0e-12)

    vessels = [
        vessel(1, "aorta_arch", {"R_poiseuille": 1.0, "C": p["fCAo"] * c("Cao"), "L": 0.0}),
        vessel(101, "upper_art", {"R_poiseuille": r("Ruba"), "C": c("Caub"), "L": l("Luba")}),
        vessel(102, "upper_cap", {"R_poiseuille": r("Rub"), "C": c("Cvub"), "L": 0.0}),
        vessel(103, "upper_vein", {"R_poiseuille": r("Rubv"), "C": 0.0, "L": 0.0}),
        vessel(104, "svc", {"R_poiseuille": r("Rsvc"), "C": c("Csvc"), "L": l("Luc")}),
        vessel(2, "aorta_th", {"R_poiseuille": r("Rthao"), "C": c("Cthao"), "L": l("Lthao")}),
        vessel(4, "aorta_ab", {"R_poiseuille": r("Rabao"), "C": c("Cabao"), "L": l("Labao")}),
        vessel(301, "leg_art", {"R_poiseuille": r("Rlega"), "C": c("Clega"), "L": l("Llega")}),
        vessel(302, "leg_cap", {"R_poiseuille": r("Rlegc"), "C": c("Clegv"), "L": 0.0}),
        vessel(303, "leg_vein", {"R_poiseuille": r("Rlegv"), "C": c("Cabivc"), "L": 0.0}),
        vessel(5, "ivc_ab", {"R_poiseuille": r("Rabivc"), "C": 0.0, "L": 0.0}),
        vessel(304, "liver_art", {"R_poiseuille": r("Rla"), "C": c("Cl"), "L": 0.0}),
        vessel(305, "liver_vein", {"R_poiseuille": r("Rlv"), "C": 0.0, "L": 0.0}),
        vessel(306, "kid_art", {"R_poiseuille": r("Rka"), "C": c("Ck"), "L": 0.0}),
        vessel(307, "kid_vein", {"R_poiseuille": r("Rkv"), "C": 0.0, "L": 0.0}),
        vessel(308, "int_art", {"R_poiseuille": r("Ria"), "C": c("Ci"), "L": 0.0}),
        vessel(309, "int_vein", {"R_poiseuille": r("Riv"), "C": 0.0, "L": 0.0}),
        vessel(310, "ivc", {"R_poiseuille": r("Rthivc"), "C": c("Cthivc"), "L": l("Llc")}),
        vessel(3, "Fontan", {"R_poiseuille": 1.0, "C": c("Ctcpc"), "L": 0.0}),
        vessel(201, "left_pul_art", {"R_poiseuille": r("Rapoll"), "C": c("Clungl"), "L": l("Lapoll")}),
        vessel(202, "left_pul_cap", {"R_poiseuille": r("Rlungl"), "C": c("Cvpl"), "L": l("Llungl")}),
        vessel(203, "left_pul_vein", {"R_poiseuille": r("Rvpl"), "C": 0.0, "L": 0.0}),
        vessel(204, "right_pul_art", {"R_poiseuille": r("Rapolr"), "C": c("Clungr"), "L": l("Lapolr")}),
        vessel(205, "right_pul_cap", {"R_poiseuille": r("Rlungr"), "C": c("Cvpr"), "L": l("Llungr")}),
        vessel(206, "right_pul_vein", {"R_poiseuille": r("Rvpr"), "C": 0.0, "L": 0.0}),
        vessel(401, "collateral_left", {"R_poiseuille": rcoll_l * P_CONV, "C": 0.0, "L": 0.0}),
        vessel(402, "collateral_right", {"R_poiseuille": rcoll_r * P_CONV, "C": 0.0, "L": 0.0}),
    ]

    junctions = [
        {"junction_name": "J3", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [1], "outlet_vessels": [101, 2]},
        {"junction_name": "J10", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [101], "outlet_vessels": [102]},
        {"junction_name": "J11", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [102], "outlet_vessels": [103]},
        {"junction_name": "J12", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [103], "outlet_vessels": [104]},
        {"junction_name": "J13", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [2], "outlet_vessels": [304, 306, 4, 401, 402]},
        {"junction_name": "J37", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [4], "outlet_vessels": [308, 301]},
        {"junction_name": "J31", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [301], "outlet_vessels": [302]},
        {"junction_name": "J32", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [302], "outlet_vessels": [303]},
        {"junction_name": "J38", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [303], "outlet_vessels": [5]},
        {"junction_name": "J35", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [308], "outlet_vessels": [309]},
        {"junction_name": "J33", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [304, 309], "outlet_vessels": [305]},
        {"junction_name": "J34", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [306], "outlet_vessels": [307]},
        {"junction_name": "J36", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [5, 305, 307], "outlet_vessels": [310]},
        {"junction_name": "J6", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [104, 310], "outlet_vessels": [3]},
        {"junction_name": "J7", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [3], "outlet_vessels": [201, 204]},
        {"junction_name": "J20", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [201], "outlet_vessels": [202]},
        {"junction_name": "J21", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [202, 401], "outlet_vessels": [203]},
        {"junction_name": "J22", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [204], "outlet_vessels": [205]},
        {"junction_name": "J23", "junction_type": "NORMAL_JUNCTION", "inlet_vessels": [205, 402], "outlet_vessels": [206]},
        {
            "junction_name": "J2",
            "junction_type": "NORMAL_JUNCTION",
            "inlet_blocks": ["left_pul_vein", "right_pul_vein"],
            "outlet_blocks": ["left_atrium"],
        },
    ]

    # Linear slopes around the initial volumes. The vascular network is exact;
    # these are the closest built-in sv0D counterparts to the genBC nonlinear
    # chamber laws.
    dv_la = X0_MMHG[26] - p["Vsa0"]
    dv_lv = X0_MMHG[28] - p["Vsv0"]
    tuned = {
        "la_emax": 1.3271785,
        "la_epass": 0.6215793,
        "lv_emax": 0.31196982,
        "lv_epass": 2.08934798,
    }
    la_epass = tuned["la_epass"] * p["csa"] * p["dsa"] * pow(2.718281828459045, p["dsa"] * dv_la) * P_CONV
    lv_epass = tuned["lv_epass"] * p["csv"] * p["dsv"] * pow(2.718281828459045, p["dsv"] * dv_lv) * P_CONV
    la_emax = tuned["la_emax"] * p["fS"] / p["Ccsa"] * P_CONV
    lv_emax = tuned["lv_emax"] * 0.5 * p["fS"] * (2.0 * p["a"] * max(dv_lv, 70.0) + p["b"]) * P_CONV

    chambers = [
        {
            "type": "LinearElastanceChamber",
            "name": "left_atrium",
            "values": {"Emax": la_emax, "Epass": la_epass, "Vrest": 0.0},
            "activation_function": {
                "type": "piecewise_cosine",
                "contract_start": 0.54336693835,
                "relax_start": 0.6407169870250001,
                "contract_duration": 0.09735004867500001,
                "relax_duration": 0.09735004867500001,
            },
        },
        {
            "type": "LinearElastanceChamber",
            "name": "left_ventricle",
            "values": {"Emax": lv_emax, "Epass": lv_epass, "Vrest": 0.0},
            "activation_function": {
                "type": "piecewise_cosine",
                "contract_start": 0.65400449366875,
                "relax_start": 0.7757101240772303,
                "contract_duration": 0.12170563040848034,
                "relax_duration": 0.06823196456026964,
            },
        },
    ]

    valve_rmin = 6.665
    valve_rmax = 6665.0
    valves = [
        {
            "type": "PiecewiseValve",
            "name": "mitral",
            "params": {
                "Rmax": valve_rmax,
                "Rmin": valve_rmin,
                "upstream_block": "left_atrium",
                "downstream_block": "left_ventricle",
            },
        },
        {
            "type": "PiecewiseValve",
            "name": "aortic",
            "params": {
                "Rmax": valve_rmax,
                "Rmin": valve_rmin,
                "upstream_block": "left_ventricle",
                "downstream_block": "aorta_arch",
            },
        },
    ]

    def sv0d_piecewise_activation(t: float, activation: dict[str, float]) -> float:
        cardiac_period = 0.666667
        phase = t % cardiac_period
        phase_unwrapped = phase + cardiac_period if phase < activation["contract_start"] else phase
        if phase_unwrapped < activation["contract_start"]:
            return 0.0
        if phase_unwrapped < activation["contract_start"] + activation["contract_duration"]:
            tau = (phase_unwrapped - activation["contract_start"]) / activation["contract_duration"]
            return 0.5 * (1.0 - math.cos(math.pi * tau))
        if phase_unwrapped < activation["relax_start"]:
            return 1.0
        if phase_unwrapped < activation["relax_start"] + activation["relax_duration"]:
            tau = (phase_unwrapped - activation["relax_start"]) / activation["relax_duration"]
            return 0.5 * (1.0 + math.cos(math.pi * tau))
        return 0.0

    def sv0d_chamber_pressure_mmhg(chamber: dict, volume: float, t: float = 0.0) -> float:
        values = chamber["values"]
        act = sv0d_piecewise_activation(t, chamber["activation_function"])
        elastance = values["Epass"] + act * (values["Emax"] - values["Epass"])
        return elastance * (volume - values["Vrest"]) / P_CONV

    def sv0d_piecewise_valve_flow_mls(p_up_mmhg: float, p_down_mmhg: float, valve: dict) -> float:
        dp_cgs = (p_up_mmhg - p_down_mmhg) * P_CONV
        resistance = valve["params"]["Rmin"] if dp_cgs >= 0.0 else valve["params"]["Rmax"]
        return dp_cgs / resistance

    p_la0 = sv0d_chamber_pressure_mmhg(chambers[0], X0_MMHG[26])
    p_lv0 = sv0d_chamber_pressure_mmhg(chambers[1], X0_MMHG[28])
    q_mitral0 = sv0d_piecewise_valve_flow_mls(p_la0, p_lv0, valves[0])
    q_aortic0 = sv0d_piecewise_valve_flow_mls(p_lv0, X0_MMHG[0], valves[1])
    q_la_in0 = (X0_MMHG[21] - p_la0) / p["Rvpr"] + (X0_MMHG[25] - p_la0) / p["Rvpl"]

    initial_condition = {
        "Vc:left_atrium": X0_MMHG[26],
        "Vc:left_ventricle": X0_MMHG[28],
        "flow:J2:left_atrium": q_la_in0,
        "flow:left_atrium:mitral": 0.0,
        "flow:mitral:left_ventricle": 0.0,
        "flow:left_ventricle:aortic": 0.0,
        "flow:aortic:aorta_arch": 0.0,
        "flow:aorta_arch:J3": X0_MMHG[1] + X0_MMHG[5],
        "flow:J3:upper_art": X0_MMHG[1],
        "flow:J3:aorta_th": X0_MMHG[5],
        "flow:upper_art:J10": X0_MMHG[1],
        "flow:J10:upper_cap": (X0_MMHG[2] - X0_MMHG[3]) / p["Rub"],
        "flow:upper_cap:J11": (X0_MMHG[2] - X0_MMHG[3]) / p["Rub"],
        "flow:J11:upper_vein": (X0_MMHG[3] - X0_MMHG[4]) / p["Rubv"],
        "flow:upper_vein:J12": (X0_MMHG[3] - X0_MMHG[4]) / p["Rubv"],
        "flow:J12:svc": abs(X0_MMHG[30]),
        "flow:svc:J6": abs(X0_MMHG[30]),
        "flow:aorta_th:J13": X0_MMHG[5],
        "flow:J13:aorta_ab": X0_MMHG[7],
        "flow:aorta_ab:J37": X0_MMHG[7],
        "flow:J37:leg_art": X0_MMHG[9],
        "flow:leg_art:J31": X0_MMHG[9],
        "flow:J31:leg_cap": (X0_MMHG[10] - X0_MMHG[11]) / p["Rlegc"],
        "flow:leg_cap:J32": (X0_MMHG[10] - X0_MMHG[11]) / p["Rlegc"],
        "flow:J32:leg_vein": (X0_MMHG[11] - X0_MMHG[12]) / p["Rlegv"],
        "flow:leg_vein:J38": (X0_MMHG[11] - X0_MMHG[12]) / p["Rlegv"],
        "flow:J38:ivc_ab": (X0_MMHG[12] - X0_MMHG[13]) / p["Rabivc"],
        "flow:ivc_ab:J36": (X0_MMHG[12] - X0_MMHG[13]) / p["Rabivc"],
        "flow:J13:liver_art": (X0_MMHG[6] - X0_MMHG[14]) / p["Rla"],
        "flow:liver_art:J33": (X0_MMHG[6] - X0_MMHG[14]) / p["Rla"],
        "flow:liver_vein:J36": (X0_MMHG[14] - X0_MMHG[13]) / p["Rlv"],
        "flow:J33:liver_vein": (X0_MMHG[14] - X0_MMHG[13]) / p["Rlv"],
        "flow:J13:kid_art": (X0_MMHG[6] - X0_MMHG[16]) / p["Rka"],
        "flow:kid_art:J34": (X0_MMHG[6] - X0_MMHG[16]) / p["Rka"],
        "flow:J34:kid_vein": (X0_MMHG[16] - X0_MMHG[13]) / p["Rkv"],
        "flow:kid_vein:J36": (X0_MMHG[16] - X0_MMHG[13]) / p["Rkv"],
        "flow:J37:int_art": (X0_MMHG[8] - X0_MMHG[15]) / p["Ria"],
        "flow:int_art:J35": (X0_MMHG[8] - X0_MMHG[15]) / p["Ria"],
        "flow:J35:int_vein": (X0_MMHG[15] - X0_MMHG[14]) / p["Riv"],
        "flow:int_vein:J33": (X0_MMHG[15] - X0_MMHG[14]) / p["Riv"],
        "flow:J36:ivc": abs(X0_MMHG[29]),
        "flow:ivc:J6": abs(X0_MMHG[29]),
        "flow:J6:Fontan": abs(X0_MMHG[29]) + abs(X0_MMHG[30]),
        "flow:Fontan:J7": X0_MMHG[18] + X0_MMHG[22],
        "flow:J7:left_pul_art": X0_MMHG[22],
        "flow:left_pul_art:J20": X0_MMHG[22],
        "flow:J20:left_pul_cap": X0_MMHG[24],
        "flow:left_pul_cap:J21": X0_MMHG[24],
        "flow:J21:left_pul_vein": (X0_MMHG[25] - 10.0) / p["Rvpl"],
        "flow:left_pul_vein:J2": (X0_MMHG[25] - 10.0) / p["Rvpl"],
        "flow:J7:right_pul_art": X0_MMHG[18],
        "flow:right_pul_art:J22": X0_MMHG[18],
        "flow:J22:right_pul_cap": X0_MMHG[20],
        "flow:right_pul_cap:J23": X0_MMHG[20],
        "flow:J23:right_pul_vein": (X0_MMHG[21] - 10.0) / p["Rvpr"],
        "flow:right_pul_vein:J2": (X0_MMHG[21] - 10.0) / p["Rvpr"],
        "flow:J13:collateral_left": (1.0 - split) * (X0_MMHG[6] - X0_MMHG[25]) / p["Rcoll"],
        "flow:collateral_left:J21": (1.0 - split) * (X0_MMHG[6] - X0_MMHG[25]) / p["Rcoll"],
        "flow:J13:collateral_right": split * (X0_MMHG[6] - X0_MMHG[21]) / p["Rcoll"],
        "flow:collateral_right:J23": split * (X0_MMHG[6] - X0_MMHG[21]) / p["Rcoll"],
    }

    pressure_ic = {
        "pressure:aorta_arch:J3": X0_MMHG[0],
        "pressure:J3:upper_art": X0_MMHG[0],
        "pressure:J10:upper_cap": X0_MMHG[2],
        "pressure:J11:upper_vein": X0_MMHG[3],
        "pressure:J12:svc": X0_MMHG[4],
        "pressure:svc:J6": X0_MMHG[17],
        "pressure:aorta_th:J13": X0_MMHG[6],
        "pressure:J13:aorta_ab": X0_MMHG[6],
        "pressure:aorta_ab:J37": X0_MMHG[8],
        "pressure:J37:leg_art": X0_MMHG[8],
        "pressure:J31:leg_cap": X0_MMHG[10],
        "pressure:J32:leg_vein": X0_MMHG[11],
        "pressure:J38:ivc_ab": X0_MMHG[12],
        "pressure:ivc_ab:J36": X0_MMHG[13],
        "pressure:J13:liver_art": X0_MMHG[6],
        "pressure:J33:liver_vein": X0_MMHG[14],
        "pressure:liver_vein:J36": X0_MMHG[13],
        "pressure:J13:kid_art": X0_MMHG[6],
        "pressure:J34:kid_vein": X0_MMHG[16],
        "pressure:kid_vein:J36": X0_MMHG[13],
        "pressure:J37:int_art": X0_MMHG[8],
        "pressure:J35:int_vein": X0_MMHG[15],
        "pressure:int_vein:J33": X0_MMHG[14],
        "pressure:J36:ivc": X0_MMHG[13],
        "pressure:ivc:J6": X0_MMHG[17],
        "pressure:J6:Fontan": X0_MMHG[17],
        "pressure:Fontan:J7": X0_MMHG[17],
        "pressure:J7:left_pul_art": X0_MMHG[17],
        "pressure:left_pul_art:J20": X0_MMHG[23],
        "pressure:J20:left_pul_cap": X0_MMHG[23],
        "pressure:left_pul_cap:J21": X0_MMHG[25],
        "pressure:J21:left_pul_vein": X0_MMHG[25],
        "pressure:left_pul_vein:J2": X0_MMHG[25],
        "pressure:J7:right_pul_art": X0_MMHG[17],
        "pressure:right_pul_art:J22": X0_MMHG[19],
        "pressure:J22:right_pul_cap": X0_MMHG[19],
        "pressure:right_pul_cap:J23": X0_MMHG[21],
        "pressure:J23:right_pul_vein": X0_MMHG[21],
        "pressure:right_pul_vein:J2": X0_MMHG[21],
        "pressure:J2:left_atrium": p_la0,
        "pressure:left_atrium:mitral": p_la0,
        "pressure:mitral:left_ventricle": p_lv0,
        "pressure:left_ventricle:aortic": p_lv0,
        "pressure:aortic:aorta_arch": X0_MMHG[0],
        "pressure:J13:collateral_left": X0_MMHG[6],
        "pressure:collateral_left:J21": X0_MMHG[25],
        "pressure:J13:collateral_right": X0_MMHG[6],
        "pressure:collateral_right:J23": X0_MMHG[21],
    }
    initial_condition.update({k: v * P_CONV for k, v in pressure_ic.items()})

    return {
        "simulation_parameters": {
            "number_of_cardiac_cycles": 20,
            "number_of_time_pts_per_cardiac_cycle": 250,
            "cardiac_period": 0.666667,
            "output_variable_based": True,
            "output_all_cycles": True,
            "steady_initial": False,
        },
        "vessels": vessels,
        "junctions": junctions,
        "boundary_conditions": [],
        "chambers": chambers,
        "valves": valves,
        "initial_condition": initial_condition,
        "_notes": {
            "source": "passive_fluid_truncate genBC closed-loop Fontan LPN",
            "unit_conversion": "R and L multiplied by 1334; C divided by 1334; flows remain mL/s.",
            "heart_translation": "This case uses the older explicit-heart svZeroDSolver schema with LinearElastanceChamber and PiecewiseValve. The vascular network remains a direct genBC translation. Chamber stiffnesses start from genBC pressure-volume slopes and are tuned with fixed vascular parameters to reduce differences in max Pao, Psa, Psv, and collateral flow. Applied multipliers are LA Emax 1.3271785, LA Epass 0.6215793, SV Emax 0.31196982 times the previously halved slope-derived Emax, and SV Epass 2.08934798. Chamber Vrest is set to 0.0. Valve Rmax is set to 6665.0 while keeping Rmin at 6.665. cardiac_period is set to 0.666667 s to match genBC HR=90. Activation times are in seconds, converted from the original genBC phase coordinates by multiplying by cardiac_period. GenBC activation wraps across t=0, so sv0D relax_start/relax_end may intentionally exceed one cardiac_period.",
            "collateral_topology": "Collateral vessels originate from the thoracic aorta branch junction J13 and merge into the pulmonary venous compliance junctions: collateral_left -> J21/LPV node and collateral_right -> J23/RPV node, matching genBC Qcoll_l into x26 and Qcoll_r into x22.",
        },
    }


def code_cell(source: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source.splitlines(keepends=True),
    }


def markdown_cell(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(keepends=True)}


NOTEBOOK_CODE = r'''
import csv
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
from IPython.display import HTML, display

CASE_DIR = Path.cwd()
if CASE_DIR.name != "sv0D_LPN":
    CASE_DIR = Path("sv0D_LPN")

JSON_PATH = CASE_DIR / "Fontan_genBC_1to1_sv0D.json"
GENBC_OUT = CASE_DIR / "genBC_closed_loop_result.csv"
SV0D_OUT = CASE_DIR / "Fontan_genBC_1to1_sv0D_result.csv"
P_CONV = 1334.0
SVZEROD_SOLVER = shutil.which("svzerodsolver")
'''.strip()


GENBC_MODEL_CODE = r'''
params = {
    "BSA": 0.85, "HR": 90.0, "tsvs": 0.1013, "tsas": 0.1947, "t1": 0.0714,
    "a": 0.137, "b": 0.914, "csv": 0.102, "dsv": 0.074, "Vsv0": 18.7802,
    "Ccsa": 6.607, "csa": 0.035, "dsa": 0.062, "Vsa0": 0.5148,
    "P0ith": 0.0, "Apith": 0.0, "Trr": 0.67, "pi": 3.141593,
    "fc": 131.58, "fr": 32.895, "Rmyo": 0.0085,
    "Lav": 1.66e-05, "Kav": 1.88e-04, "Kao": 1.97e-05,
    "Cao": 0.6, "Ruba": 0.761, "Luba": 6.06e-04, "Caub": 0.236,
    "Rub": 3.193, "Cvub": 3.1472, "Rubv": 0.057, "Csvc": 1.0277, "Rsvc": 0.0076,
    "Rthao": 0.047, "Lthao": 0.0018, "Cthao": 0.1773,
    "Rabao": 0.646, "Labao": 0.0018, "Cabao": 0.3,
    "Rlega": 0.929, "Llega": 0.0018, "Clega": 1.5829,
    "Rlegc": 2.36, "Clegv": 7.2032, "Rlegv": 0.548,
    "Cabivc": 6.502, "Rabivc": 0.0031, "V0abivc": 11.1282, "R0abivc": 0.04,
    "P0iab": 0.0, "Apiab": 0.0, "Cthivc": 0.72712, "Rthivc": 0.003,
    "V0thivc": 7.2898, "R0thivc": 0.01, "Llc": 0.00516, "Luc": 0.00516,
    "Klc": 0.000516, "Kuc": 0.000516, "Rla": 13.257, "Cl": 4.42, "Rlv": 0.045,
    "Rka": 9.201, "Ck": 1.83, "Rkv": 0.878, "Ria": 18.599, "Ci": 0.99, "Riv": 0.143,
    "Ctcpc": 0.5, "Rapoll": 0.073, "Rapolr": 0.02057, "Lapoll": 0.00227, "Lapolr": 0.00227,
    "Clungl": 0.8335, "Clungr": 1.0043, "Rlungl": 0.121, "Rlungr": 0.05569,
    "Llungl": 0.00153, "Llungr": 0.00153, "Cvpl": 1.5904, "Cvpr": 1.8495,
    "Rvpl": 0.034, "Rvpr": 0.01439, "Clpa": 1.3111, "Crpa": 1.3111,
    "coll": 1.0, "Rcoll": 5.636, "coll_split": 0.5, "Cub": 4.184140958,
    "fS": 0.797, "fD": 1.0, "fCAo": 1.21,
}

x0 = np.array([
    84.91167, 39.37003, 55.49430, 12.3695, 12.39335, 27.36710, 74.04834, 24.2435,
    55.40978, 10.50352, 45.0045, 18.56190, 12.61708, 12.59937, 12.06896, 11.8046,
    16.02366, 12.25075, 11.809034, 11.05807, 10.000172, 9.849035, 10.724361,
    11.49486, 14.955302, 10.07034, 72.11296, 0.0, 43.5820, -26.0, -15.0, 0.0
], dtype=float)

elastance_coeffs = np.array([
    0.00, 0.000, 0.0000, 0.0000, 11.56, -0.190, 0.0000, 0.0000,
    0.30, 4.140, 0.4938, 0.0196, -25.42, 4.270, 1.5450, 0.1466,
    16.86, -5.275, 1.4200, 0.3568, -3.29, 1.045, 0.8912, 0.4848,
    1.95, -0.200, 0.9969, 0.6061, -4.31, 0.545, 1.0400, 0.7314,
    -0.11, -1.070, 0.9738, 0.8615, -30.56, -1.095, 0.7031, 0.9663,
    23.76, -12.560, -1.0037, 0.9774, 32.02, -3.630, -3.0281, 0.7021,
    -15.49, 8.365, -2.4356, 0.3294, 0.35, 2.575, -1.0681, 0.1254,
    -7.67, 2.700, -0.4088, 0.0328, 3.04, -0.190, -0.0950, 0.0089,
    0.00, 0.000, 0.0000, 0.0000,
])

def elastance(t, p=params, coeffs=elastance_coeffs):
    tc = 60.0 / p["HR"]
    tmv = ((t + 0.125 * p["tsvs"]) % tc) / p["tsvs"]
    if tmv < 1.875:
        seg = int(np.floor(tmv / 0.125))
        dt = tmv - seg * 0.125
        i = 4 * seg
        return ((coeffs[i] * dt + coeffs[i + 1]) * dt + coeffs[i + 2]) * dt + coeffs[i + 3]
    return 0.0

def rhs(x, p=params):
    f = np.zeros(32)
    y = {}
    t = x[31]
    tc = 60.0 / p["HR"]
    tr = 4.0 * tc

    tma = (t + p["tsas"] - p["t1"]) % tc
    aa = (1.0 - np.cos(2.0 * p["pi"] * tma / p["tsas"])) / 2.0 if tma < p["tsas"] else 0.0
    dvsa = x[26] - p["Vsa0"]
    psa = p["fS"] * aa * dvsa / p["Ccsa"] + p["fD"] * p["csa"] * (np.exp(p["dsa"] * dvsa) - 1.0)

    dvsv = x[28] - p["Vsv0"]
    psvp = p["csv"] * (np.exp(p["dsv"] * dvsv) - 1.0)
    psva = (p["a"] * dvsv + p["b"]) * dvsv
    psv = p["fS"] * elastance(t, p) * psva + p["fD"] * psvp

    tmr = t % tr
    tsr = p["Trr"] * tr
    if tmr < tsr:
        dar = (p["pi"] / tsr) * np.sin(2.0 * p["pi"] * tmr / tsr)
    else:
        dar = 0.0
    dpith = p["Apith"] * dar

    if x[27] < 0.0:
        x = x.copy()
        x[27] = 0.0

    qao = (np.sqrt(p["Rmyo"] ** 2 + 4.0 * p["Kao"] * (psv - x[0])) - p["Rmyo"]) / (2.0 * p["Kao"]) if psv > x[0] else 0.0
    qub = (x[2] - x[3]) / p["Rub"]
    qvub = (x[3] - x[4]) / p["Rubv"]
    qsvc = x[30]
    qlegc = (x[10] - x[11]) / p["Rlegc"]
    qlegv = (abs(x[11] - x[12]) + x[11] - x[12]) / (2.0 * p["Rlegv"])
    qabivc = (x[12] - x[13]) / p["Rabivc"]
    qthivc = x[29]
    qla = (x[6] - x[14]) / p["Rla"]
    qlv = (x[14] - x[13]) / p["Rlv"]
    qia = (x[8] - x[15]) / p["Ria"]
    qiv = (x[15] - x[14]) / p["Riv"]
    qka = (x[6] - x[16]) / p["Rka"]
    qkv = (x[16] - x[13]) / p["Rkv"]
    qvpr = (x[21] - psa) / p["Rvpr"]
    qvpl = (x[25] - psa) / p["Rvpl"]
    qcoll_r = p["coll_split"] * (x[6] - x[21]) / p["Rcoll"]
    qcoll_l = (1.0 - p["coll_split"]) * (x[6] - x[25]) / p["Rcoll"]
    qcoll = qcoll_r + qcoll_l

    f[0] = (qao - x[1] - x[5]) / (p["fCAo"] * p["Cao"]) + dpith
    f[1] = (x[0] - x[1] * p["Ruba"] - x[2]) / p["Luba"]
    f[2] = (x[1] - qub) / p["Caub"]
    f[3] = (qub - qvub) / p["Cvub"]
    f[4] = (qvub - qsvc) / p["Csvc"] + dpith
    f[5] = (x[0] - x[5] * p["Rthao"] - x[6]) / p["Lthao"]
    f[6] = (x[5] - qcoll - x[7] - qla - qka) / p["Cthao"] + dpith
    f[7] = (x[6] - x[7] * p["Rabao"] - x[8]) / p["Labao"]
    f[8] = (x[7] - qia - x[9]) / p["Cabao"]
    f[9] = (x[8] - x[9] * p["Rlega"] - x[10]) / p["Llega"]
    f[10] = (x[9] - qlegc) / p["Clega"]
    f[11] = (qlegc - qlegv) / p["Clegv"]
    f[12] = (qlegv - qabivc) / p["Cabivc"]
    f[13] = (qabivc + qkv + qlv - qthivc) / p["Cthivc"] + dpith
    f[14] = (qiv + qla - qlv) / p["Cl"]
    f[15] = (qia - qiv) / p["Ci"]
    f[16] = (qka - qkv) / p["Ck"]
    f[17] = (qsvc + qthivc - x[18] - x[22]) / p["Ctcpc"] + dpith
    f[18] = (x[17] - x[18] * p["Rapolr"] - x[19]) / p["Lapolr"]
    f[19] = (x[18] - x[20]) / p["Clungr"] + dpith
    f[20] = (x[19] - x[20] * p["Rlungr"] - x[21]) / p["Llungr"]
    f[21] = (x[20] + qcoll_r - qvpr) / p["Cvpr"] + dpith
    f[22] = (x[17] - x[22] * p["Rapoll"] - x[23]) / p["Lapoll"]
    f[23] = (x[22] - x[24]) / p["Clungl"] + dpith
    f[24] = (x[23] - x[24] * p["Rlungl"] - x[25]) / p["Llungl"]
    f[25] = (x[24] + qcoll_l - qvpl) / p["Cvpl"] + dpith
    f[26] = qvpr + qvpl - x[27]
    f[27] = 0.0 if (x[27] == 0.0 and psa < (psv - p["Rmyo"] * qao)) else (psa - psv - p["Kav"] * x[27] ** 2) / p["Lav"]
    f[28] = x[27] - qao
    f[29] = (x[13] - x[29] * p["Rthivc"] - x[17]) / p["Llc"]
    f[30] = (x[4] - x[30] * p["Rsvc"] - x[17]) / p["Luc"]
    f[31] = 1.0

    y.update(Psa=psa, Psv=psv, Qao=qao, Qcoll_l=qcoll_l, Qcoll_r=qcoll_r, Qcoll=qcoll)
    return f, y

def run_genbc_reference():
    steps = 250
    substeps = 40
    ncyc = 20
    dt = (60.0 / params["HR"]) / (steps * substeps)
    total_steps = ncyc * steps * substeps
    out = []
    x = x0.copy()

    for i in range(total_steps + 1):
        if i % substeps == 0:
            _, y = rhs(x)
            out.append([
                x[31], *x[:31], y["Psa"], y["Psv"], y["Qao"],
                y["Qcoll_l"], y["Qcoll_r"], y["Qcoll"]
            ])
        if i == total_steps:
            break
        k1, _ = rhs(x)
        k2, _ = rhs(x + 0.5 * dt * k1)
        k3, _ = rhs(x + 0.5 * dt * k2)
        k4, _ = rhs(x + dt * k3)
        x = x + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

    cols = ["time"] + [f"x{i}" for i in range(1, 32)] + ["Psa", "Psv", "Qao", "Qcoll_l", "Qcoll_r", "Qcoll"]
    with GENBC_OUT.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(cols)
        writer.writerows(out)
    arr = np.asarray(out, dtype=float)
    return {name: arr[:, i] for i, name in enumerate(cols)}
'''.strip()


RUN_AND_COMPARE_CODE = r'''
genbc = run_genbc_reference()
if SVZEROD_SOLVER is None:
    raise FileNotFoundError("svzerodsolver was not found on PATH")
subprocess.run([SVZEROD_SOLVER, str(JSON_PATH), str(SV0D_OUT)], check=True)

sv0d = {}
with SV0D_OUT.open(newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        entry = sv0d.setdefault(row["name"], {"time": [], "y": []})
        entry["time"].append(float(row["time"]))
        entry["y"].append(float(row["y"]))
for entry in sv0d.values():
    entry["time"] = np.asarray(entry["time"])
    entry["y"] = np.asarray(entry["y"])

def pick(names):
    for name in names:
        if name in sv0d:
            return sv0d[name]
    raise KeyError(f"None of these columns are present: {names}")

with JSON_PATH.open() as f:
    sv0d_config_for_time = json.load(f)
sv0d_period = sv0d_config_for_time["simulation_parameters"].get("cardiac_period")

compare = {
    "time": genbc["time"],
    "genBC_Pao": genbc["x1"],
    "genBC_Pfontan": genbc["x18"],
    "genBC_Qlpa": genbc["x23"],
    "genBC_Qrpa": genbc["x19"],
    "genBC_Qivc": genbc["x30"],
    "genBC_Qsvc": genbc["x31"],
    "genBC_Qao": genbc["Qao"],
    "genBC_Qmitral": genbc["x28"],
    "genBC_Pla": genbc["Psa"],
    "genBC_Plv": genbc["Psv"],
    "genBC_Vla": genbc["x27"],
    "genBC_Vlv": genbc["x29"],
    "genBC_Qcollateral": genbc["Qcoll"],
}

def interp(entry, scale=1.0):
    if sv0d_period is None:
        sv0d_physical_time = entry["time"] * (60.0 / params["HR"])
    else:
        sv0d_physical_time = entry["time"]
    return np.interp(compare["time"], sv0d_physical_time, entry["y"] * scale)

compare["sv0D_Pao"] = interp(pick(["pressure:aorta_arch:J3", "pressure:aortic:aorta_arch"]), 1.0 / P_CONV)
compare["sv0D_Pfontan"] = interp(pick(["pressure:Fontan:J7", "pressure:J6:Fontan"]), 1.0 / P_CONV)
compare["sv0D_Qlpa"] = interp(pick(["flow:J7:left_pul_art", "flow:Fontan:J7"]))
compare["sv0D_Qrpa"] = interp(pick(["flow:J7:right_pul_art"]))
compare["sv0D_Qivc"] = interp(pick(["flow:ivc:J6", "flow:J36:ivc"]))
compare["sv0D_Qsvc"] = interp(pick(["flow:svc:J6", "flow:J12:svc"]))
compare["sv0D_Qao"] = interp(pick(["flow:aortic:aorta_arch", "flow:left_ventricle:aortic"]))
compare["sv0D_Qmitral"] = interp(pick(["flow:left_atrium:mitral", "flow:mitral:left_ventricle"]))
compare["sv0D_Vla"] = interp(pick(["Vc:left_atrium"]))
compare["sv0D_Vlv"] = interp(pick(["Vc:left_ventricle"]))
compare["sv0D_Pla"] = interp(pick(["pressure:J2:left_atrium"]), 1.0 / P_CONV)
compare["sv0D_Plv"] = interp(pick(["pressure:left_ventricle:aortic", "pressure:mitral:left_ventricle"]), 1.0 / P_CONV)
compare["sv0D_Qcollateral"] = (
    interp(pick(["flow:collateral_left:J21", "flow:J13:collateral_left"]))
    + interp(pick(["flow:collateral_right:J23", "flow:J13:collateral_right"]))
)

cols = list(compare)
with (CASE_DIR / "comparison_timeseries.csv").open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(cols)
    for i in range(len(compare["time"])):
        writer.writerow([compare[name][i] for name in cols])

display(HTML(f"<p>Wrote <code>{GENBC_OUT}</code>, <code>{SV0D_OUT}</code>, and <code>{CASE_DIR / 'comparison_timeseries.csv'}</code>.</p>"))
'''.strip()


METRICS_AND_PLOTS_CODE = r'''
tc = 60.0 / params["HR"]
last = compare["time"] >= (compare["time"].max() - tc)

def mean_col(name):
    return float(np.mean(compare[name][last]))

def max_col(name):
    return float(np.max(compare[name][last]))

def l2_rel(a, b):
    aa = compare[a][last]
    bb = compare[b][last]
    return np.linalg.norm(aa - bb) / (np.linalg.norm(bb) + 1.0e-12)

metrics = [
    ["Pressure", "Pao", "Mean", "mmHg", mean_col("genBC_Pao"), mean_col("sv0D_Pao")],
    ["Pressure", "Pao", "Max", "mmHg", max_col("genBC_Pao"), max_col("sv0D_Pao")],
    ["Pressure", "Fontan", "Mean", "mmHg", mean_col("genBC_Pfontan"), mean_col("sv0D_Pfontan")],
    ["Pressure", "Psa / Pla", "Mean", "mmHg", mean_col("genBC_Pla"), mean_col("sv0D_Pla")],
    ["Pressure", "Psa / Pla", "Max", "mmHg", max_col("genBC_Pla"), max_col("sv0D_Pla")],
    ["Pressure", "Psv / Plv", "Mean", "mmHg", mean_col("genBC_Plv"), mean_col("sv0D_Plv")],
    ["Pressure", "Psv / Plv", "Max", "mmHg", max_col("genBC_Plv"), max_col("sv0D_Plv")],
    ["Flow", "Qao", "Mean", "mL/s", mean_col("genBC_Qao"), mean_col("sv0D_Qao")],
    ["Flow", "Qao", "Max", "mL/s", max_col("genBC_Qao"), max_col("sv0D_Qao")],
    ["Flow", "Collateral", "Mean", "mL/s", mean_col("genBC_Qcollateral"), mean_col("sv0D_Qcollateral")],
    ["Flow", "Collateral", "Max", "mL/s", max_col("genBC_Qcollateral"), max_col("sv0D_Qcollateral")],
    ["Flow", "LPA", "Mean", "mL/s", mean_col("genBC_Qlpa"), mean_col("sv0D_Qlpa")],
    ["Flow", "RPA", "Mean", "mL/s", mean_col("genBC_Qrpa"), mean_col("sv0D_Qrpa")],
    ["Flow", "IVC", "Mean", "mL/s", mean_col("genBC_Qivc"), mean_col("sv0D_Qivc")],
    ["Flow", "SVC", "Mean", "mL/s", mean_col("genBC_Qsvc"), mean_col("sv0D_Qsvc")],
]

def pct_change(s, g):
    return 100.0 * (s - g) / (abs(g) + 1.0e-12)

def fmt_value(value):
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 10:
        return f"{value:.2f}"
    return f"{value:.3f}"

def fmt_pct(value):
    return f"{100.0 * value:.1f}%"

def error_class(err):
    if err < 0.10:
        return "good"
    if err < 0.25:
        return "moderate"
    return "large"

table_rows = []
for group, variable, statistic, unit, genbc_value, sv0d_value in metrics:
    relative_error = abs(sv0d_value - genbc_value) / (abs(genbc_value) + 1.0e-12)
    table_rows.append({
        "group": group,
        "metric": f"{variable} {statistic.lower()}",
        "unit": unit,
        "genBC": genbc_value,
        "sv0D": sv0d_value,
        "relative_error": relative_error,
    })

html_rows = "\n".join(
    "<tr>"
    f"<td class='group {row['group'].lower()}'>{row['group']}</td>"
    f"<td class='metric'>{row['metric']}</td>"
    f"<td class='unit'>{row['unit']}</td>"
    f"<td>{fmt_value(row['genBC'])}</td>"
    f"<td>{fmt_value(row['sv0D'])}</td>"
    f"<td><span class='badge {error_class(row['relative_error'])}'>{fmt_pct(row['relative_error'])}</span></td>"
    "</tr>"
    for row in table_rows
)

slide_table_html = f"""
<style>
.slide-table-card {{
  width: 880px;
  background: #ffffff;
  border: 1px solid #d8dee8;
  border-radius: 8px;
  overflow: hidden;
  font-family: Arial, Helvetica, sans-serif;
  color: #172033;
}}
.slide-table-title {{
  padding: 14px 18px 8px 18px;
  font-size: 20px;
  font-weight: 700;
}}
.slide-table-subtitle {{
  padding: 0 18px 12px 18px;
  color: #5b6576;
  font-size: 12px;
}}
.slide-table-card table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}}
.slide-table-card th {{
  background: #eef2f7;
  color: #334155;
  padding: 8px 10px;
  border-top: 1px solid #d8dee8;
  border-bottom: 1px solid #d8dee8;
  text-align: right;
  font-weight: 700;
}}
.slide-table-card th:nth-child(1),
.slide-table-card th:nth-child(2),
.slide-table-card th:nth-child(3) {{
  text-align: left;
}}
.slide-table-card td {{
  padding: 7px 10px;
  border-bottom: 1px solid #edf1f6;
  text-align: right;
  white-space: nowrap;
}}
.slide-table-card tr:nth-child(even) td {{
  background: #fafbfd;
}}
.slide-table-card td.metric,
.slide-table-card td.unit,
.slide-table-card td.group {{
  text-align: left;
}}
.slide-table-card td.unit {{
  color: #64748b;
}}
.slide-table-card td.group {{
  font-weight: 700;
  border-left: 5px solid #2b6cb0;
}}
.slide-table-card td.group.flow {{
  border-left-color: #0f8b6c;
}}
.slide-table-card .badge {{
  display: inline-block;
  min-width: 52px;
  padding: 3px 7px;
  border-radius: 999px;
  text-align: center;
  font-weight: 700;
}}
.slide-table-card .badge.good {{
  background: #e7f6ef;
  color: #087049;
}}
.slide-table-card .badge.moderate {{
  background: #fff3d8;
  color: #9a5a00;
}}
.slide-table-card .badge.large {{
  background: #fde8e8;
  color: #b42318;
}}
</style>
<div class="slide-table-card">
  <div class="slide-table-title">Full 0D Last-Cycle Comparison</div>
  <div class="slide-table-subtitle">genBC baseline vs translated sv0D model. Error is the scalar relative error of each displayed mean or max value.</div>
  <table>
    <thead>
      <tr><th>Group</th><th>Metric</th><th>Unit</th><th>genBC</th><th>sv0D</th><th>Rel. error</th></tr>
    </thead>
    <tbody>{html_rows}</tbody>
  </table>
</div>
"""
display(HTML(slide_table_html))

def svg_escape(text):
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def svg_text(x, y, text, size=18, weight="400", anchor="start", fill="#172033"):
    return f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{svg_escape(text)}</text>'

svg_width = 1120
row_h = 34
top = 94
left = 28
columns = [left, 130, 385, 500, 655, 810, 1060]
svg_height = top + row_h * (len(table_rows) + 1) + 34
parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_width}" height="{svg_height}" viewBox="0 0 {svg_width} {svg_height}">',
    '<rect width="100%" height="100%" fill="#ffffff"/>',
    '<rect x="14" y="14" width="1092" height="{0}" rx="8" fill="#ffffff" stroke="#d8dee8"/>'.format(svg_height - 28),
    svg_text(34, 48, "Full 0D Last-Cycle Comparison", 24, "700"),
    svg_text(34, 72, "genBC baseline vs translated sv0D model; error = scalar relative error for each mean or max value", 13, "400", fill="#5b6576"),
    f'<rect x="{left}" y="{top}" width="{columns[-1] - left}" height="{row_h}" fill="#eef2f7"/>',
]
headers = ["Group", "Metric", "Unit", "genBC", "sv0D", "Rel. error"]
for i, header in enumerate(headers):
    anchor = "start" if i < 3 else "end"
    x = columns[i] + 10 if i < 3 else columns[i + 1] - 10
    parts.append(svg_text(x, top + 23, header, 13, "700", anchor, "#334155"))

for idx, row in enumerate(table_rows):
    y = top + row_h * (idx + 1)
    fill = "#fafbfd" if idx % 2 else "#ffffff"
    accent = "#2b6cb0" if row["group"] == "Pressure" else "#0f8b6c"
    badge_fill, badge_text = {
        "good": ("#e7f6ef", "#087049"),
        "moderate": ("#fff3d8", "#9a5a00"),
        "large": ("#fde8e8", "#b42318"),
    }[error_class(row["relative_error"])]
    parts.extend([
        f'<rect x="{left}" y="{y}" width="{columns[-1] - left}" height="{row_h}" fill="{fill}"/>',
        f'<rect x="{left}" y="{y}" width="5" height="{row_h}" fill="{accent}"/>',
        svg_text(columns[0] + 10, y + 22, row["group"], 13, "700"),
        svg_text(columns[1] + 10, y + 22, row["metric"], 13),
        svg_text(columns[2] + 10, y + 22, row["unit"], 13, fill="#64748b"),
        svg_text(columns[4] - 10, y + 22, fmt_value(row["genBC"]), 13, anchor="end"),
        svg_text(columns[5] - 10, y + 22, fmt_value(row["sv0D"]), 13, anchor="end"),
        f'<rect x="{columns[5] + 86}" y="{y + 7}" width="76" height="20" rx="10" fill="{badge_fill}"/>',
        svg_text(columns[6] - 20, y + 22, fmt_pct(row["relative_error"]), 13, "700", "end", badge_text),
    ])

parts.append("</svg>")
slide_table_svg = "\n".join(parts)
slide_table_svg_path = CASE_DIR / "fontan_0d_summary_slide_table.svg"
slide_table_svg_path.write_text(slide_table_svg)
display(HTML(f"<p>Wrote slide-ready SVG table: <code>{slide_table_svg_path}</code></p>"))

plot_window = compare["time"] >= (compare["time"].max() - 2 * tc)

def line_points(x, y, xmin, xmax, ymin, ymax, width, height, pad):
    xs = pad + (x - xmin) / max(xmax - xmin, 1.0e-12) * (width - 2 * pad)
    ys = height - pad - (y - ymin) / max(ymax - ymin, 1.0e-12) * (height - 2 * pad)
    return " ".join(f"{xx:0.1f},{yy:0.1f}" for xx, yy in zip(xs, ys))

def svg_timeseries(title, ylabel, g, s, width=420, height=250):
    pad = 42
    x = compare["time"][plot_window]
    yg = compare[g][plot_window]
    ys = compare[s][plot_window]
    ymin = float(min(np.min(yg), np.min(ys)))
    ymax = float(max(np.max(yg), np.max(ys)))
    margin = 0.08 * max(ymax - ymin, 1.0)
    ymin -= margin
    ymax += margin
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    p1 = line_points(x, yg, xmin, xmax, ymin, ymax, width, height, pad)
    p2 = line_points(x, ys, xmin, xmax, ymin, ymax, width, height, pad)
    return f"""
    <svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" style="max-width:100%;border:1px solid #ddd;background:white">
      <text x="{pad}" y="22" font-size="15" font-family="sans-serif">{title}</text>
      <line x1="{pad}" y1="{height-pad}" x2="{width-pad}" y2="{height-pad}" stroke="#555"/>
      <line x1="{pad}" y1="{pad}" x2="{pad}" y2="{height-pad}" stroke="#555"/>
      <text x="{width-pad-80}" y="{height-12}" font-size="11" font-family="sans-serif">time [s]</text>
      <text x="8" y="{pad-10}" font-size="11" font-family="sans-serif">{ylabel}</text>
      <polyline points="{p1}" fill="none" stroke="#1f77b4" stroke-width="2"/>
      <polyline points="{p2}" fill="none" stroke="#d62728" stroke-width="2" stroke-dasharray="6 4"/>
      <text x="{pad+8}" y="{height-pad-12}" font-size="11" font-family="sans-serif" fill="#1f77b4">genBC Python</text>
      <text x="{pad+108}" y="{height-pad-12}" font-size="11" font-family="sans-serif" fill="#d62728">sv0D</text>
    </svg>
    """

def svg_xy(title, xlabel, ylabel, xg, yg, xs, ys, width=420, height=250):
    pad = 42
    x1 = compare[xg][plot_window]
    y1 = compare[yg][plot_window]
    x2 = compare[xs][plot_window]
    y2 = compare[ys][plot_window]
    xmin = float(min(np.min(x1), np.min(x2)))
    xmax = float(max(np.max(x1), np.max(x2)))
    ymin = float(min(np.min(y1), np.min(y2)))
    ymax = float(max(np.max(y1), np.max(y2)))
    xmin -= 0.08 * max(xmax - xmin, 1.0)
    xmax += 0.08 * max(xmax - xmin, 1.0)
    ymin -= 0.08 * max(ymax - ymin, 1.0)
    ymax += 0.08 * max(ymax - ymin, 1.0)
    p1 = line_points(x1, y1, xmin, xmax, ymin, ymax, width, height, pad)
    p2 = line_points(x2, y2, xmin, xmax, ymin, ymax, width, height, pad)
    return f"""
    <svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" style="max-width:100%;border:1px solid #ddd;background:white">
      <text x="{pad}" y="22" font-size="15" font-family="sans-serif">{title}</text>
      <line x1="{pad}" y1="{height-pad}" x2="{width-pad}" y2="{height-pad}" stroke="#555"/>
      <line x1="{pad}" y1="{pad}" x2="{pad}" y2="{height-pad}" stroke="#555"/>
      <text x="{width-pad-90}" y="{height-12}" font-size="11" font-family="sans-serif">{xlabel}</text>
      <text x="8" y="{pad-10}" font-size="11" font-family="sans-serif">{ylabel}</text>
      <polyline points="{p1}" fill="none" stroke="#1f77b4" stroke-width="2"/>
      <polyline points="{p2}" fill="none" stroke="#d62728" stroke-width="2" stroke-dasharray="6 4"/>
    </svg>
    """

pairs = [
    ("Pao", "Pressure [mmHg]", "genBC_Pao", "sv0D_Pao"),
    ("Fontan pressure", "Pressure [mmHg]", "genBC_Pfontan", "sv0D_Pfontan"),
    ("Psa", "Pressure [mmHg]", "genBC_Pla", "sv0D_Pla"),
    ("Psv", "Pressure [mmHg]", "genBC_Plv", "sv0D_Plv"),
    ("Qao flow", "Flow [mL/s]", "genBC_Qao", "sv0D_Qao"),
    ("Collateral flow", "Flow [mL/s]", "genBC_Qcollateral", "sv0D_Qcollateral"),
    ("LPA flow", "Flow [mL/s]", "genBC_Qlpa", "sv0D_Qlpa"),
    ("RPA flow", "Flow [mL/s]", "genBC_Qrpa", "sv0D_Qrpa"),
    ("IVC flow", "Flow [mL/s]", "genBC_Qivc", "sv0D_Qivc"),
    ("SVC flow", "Flow [mL/s]", "genBC_Qsvc", "sv0D_Qsvc"),
]
plots = [svg_timeseries(*pair) for pair in pairs]
plots.append(svg_xy("Left atrium PV loop", "Volume [mL]", "Pressure [mmHg]", "genBC_Vla", "genBC_Pla", "sv0D_Vla", "sv0D_Pla"))
plots.append(svg_xy("Single ventricle PV loop", "Volume [mL]", "Pressure [mmHg]", "genBC_Vlv", "genBC_Plv", "sv0D_Vlv", "sv0D_Plv"))
display(HTML("<div style='display:grid;grid-template-columns:repeat(2,minmax(320px,1fr));gap:14px'>" + "".join(plots) + "</div>"))
'''.strip()


VALVE_ANALYSIS_CODE = r'''
tc = 60.0 / params["HR"]
last = compare["time"] >= (compare["time"].max() - tc)
plot_window = compare["time"] >= (compare["time"].max() - 2 * tc)

def trapz(y, x):
    return np.trapezoid(y, x) if hasattr(np, "trapezoid") else np.trapz(y, x)

def valve_flow_stats(label, q_name):
    t = compare["time"][last]
    q = compare[q_name][last]
    q_forward = np.maximum(q, 0.0)
    q_back = np.maximum(-q, 0.0)
    forward_volume = float(trapz(q_forward, t))
    backflow_volume = float(trapz(q_back, t))
    net_volume = float(trapz(q, t))
    backflow_fraction_forward = 100.0 * backflow_volume / (forward_volume + 1.0e-12)
    backflow_fraction_net = 100.0 * backflow_volume / (abs(net_volume) + 1.0e-12)
    negative_duration = 100.0 * float(np.mean(q < 0.0))
    return [
        label,
        float(np.mean(q)),
        float(np.max(q)),
        float(np.min(q)),
        forward_volume,
        backflow_volume,
        net_volume,
        backflow_fraction_forward,
        backflow_fraction_net,
        negative_duration,
    ]

valve_rows = [
    valve_flow_stats("genBC Qao", "genBC_Qao"),
    valve_flow_stats("sv0D aortic valve", "sv0D_Qao"),
    valve_flow_stats("genBC mitral valve", "genBC_Qmitral"),
    valve_flow_stats("sv0D mitral valve", "sv0D_Qmitral"),
]

rows = "\n".join(
    "<tr>"
    f"<td>{name}</td>"
    f"<td>{mean_q:0.4g}</td>"
    f"<td>{max_q:0.4g}</td>"
    f"<td>{min_q:0.4g}</td>"
    f"<td>{forward_vol:0.4g}</td>"
    f"<td>{back_vol:0.4g}</td>"
    f"<td>{net_vol:0.4g}</td>"
    f"<td>{back_forward:0.4g}</td>"
    f"<td>{back_net:0.4g}</td>"
    f"<td>{neg_dur:0.4g}</td>"
    "</tr>"
    for (
        name,
        mean_q,
        max_q,
        min_q,
        forward_vol,
        back_vol,
        net_vol,
        back_forward,
        back_net,
        neg_dur,
    ) in valve_rows
)
display(HTML(
    "<h3>Valve flow and backflow, last cycle</h3>"
    "<table>"
    "<thead><tr>"
    "<th>signal</th><th>mean Q [mL/s]</th><th>max Q [mL/s]</th><th>min Q [mL/s]</th>"
    "<th>forward volume [mL/beat]</th><th>backflow volume [mL/beat]</th><th>net volume [mL/beat]</th>"
    "<th>backflow / forward [%]</th><th>backflow / net [%]</th><th>time Q&lt;0 [%]</th>"
    "</tr></thead>"
    f"<tbody>{rows}</tbody></table>"
))

def valve_line_points(x, y, xmin, xmax, ymin, ymax, width, height, pad):
    xs = pad + (x - xmin) / max(xmax - xmin, 1.0e-12) * (width - 2 * pad)
    ys = height - pad - (y - ymin) / max(ymax - ymin, 1.0e-12) * (height - 2 * pad)
    return " ".join(f"{xx:0.1f},{yy:0.1f}" for xx, yy in zip(xs, ys))

def svg_valve_flow(title, g_name, s_name, genbc_label, sv0d_label, width=760, height=300):
    pad = 46
    x = compare["time"][plot_window]
    yg = compare[g_name][plot_window]
    ys = compare[s_name][plot_window]
    ymin = float(min(np.min(yg), np.min(ys), 0.0))
    ymax = float(max(np.max(yg), np.max(ys), 0.0))
    margin = 0.08 * max(ymax - ymin, 1.0)
    ymin -= margin
    ymax += margin
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    p1 = valve_line_points(x, yg, xmin, xmax, ymin, ymax, width, height, pad)
    p2 = valve_line_points(x, ys, xmin, xmax, ymin, ymax, width, height, pad)
    yzero = height - pad - (0.0 - ymin) / max(ymax - ymin, 1.0e-12) * (height - 2 * pad)
    return f"""
    <svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" style="max-width:100%;border:1px solid #ddd;background:white">
      <text x="{pad}" y="24" font-size="16" font-family="sans-serif">{title}, last two cycles</text>
      <line x1="{pad}" y1="{height-pad}" x2="{width-pad}" y2="{height-pad}" stroke="#555"/>
      <line x1="{pad}" y1="{pad}" x2="{pad}" y2="{height-pad}" stroke="#555"/>
      <line x1="{pad}" y1="{yzero:0.1f}" x2="{width-pad}" y2="{yzero:0.1f}" stroke="#444" stroke-width="1" stroke-dasharray="3 3"/>
      <text x="{width-pad-78}" y="{height-12}" font-size="11" font-family="sans-serif">time [s]</text>
      <text x="8" y="{pad-12}" font-size="11" font-family="sans-serif">Qao [mL/s]</text>
      <text x="{pad+4}" y="{yzero-5:0.1f}" font-size="10" font-family="sans-serif" fill="#444">Q=0</text>
      <polyline points="{p1}" fill="none" stroke="#1f77b4" stroke-width="2"/>
      <polyline points="{p2}" fill="none" stroke="#d62728" stroke-width="2" stroke-dasharray="6 4"/>
      <text x="{pad+8}" y="{height-pad-12}" font-size="11" font-family="sans-serif" fill="#1f77b4">{genbc_label}</text>
      <text x="{pad+122}" y="{height-pad-12}" font-size="11" font-family="sans-serif" fill="#d62728">{sv0d_label}</text>
    </svg>
    """

valve_plots = [
    svg_valve_flow("Aortic valve flow", "genBC_Qao", "sv0D_Qao", "genBC Qao", "sv0D aortic valve"),
    svg_valve_flow("Mitral valve flow", "genBC_Qmitral", "sv0D_Qmitral", "genBC mitral", "sv0D mitral valve"),
]
display(HTML("<div style='display:grid;grid-template-columns:1fr;gap:14px'>" + "".join(valve_plots) + "</div>"))
'''.strip()


PARAMETER_IC_AUDIT_CODE = r'''
import re

source_root = CASE_DIR.parent
param_file = source_root / "passive_fluid_truncate/genBC_QP_Baretta/include/parameters_final.f"
init_file = source_root / "passive_fluid_truncate/genBC_QP_Baretta/include/initial_values_final.f"
sv0d_config = json.loads(JSON_PATH.read_text())
sv0d_ic = sv0d_config["initial_condition"]

def parse_fortran_assignments(path):
    text = path.read_text()
    found = {}
    pattern = r"([A-Za-z][A-Za-z0-9_]*)\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?)"
    for name, value in re.findall(pattern, text):
        found[name.lower()] = float(value.replace("D", "E").replace("d", "e"))
    return found

def parse_fortran_initial_values(path):
    text = path.read_text().split("(/", 1)[1].split("/)", 1)[0]
    return np.array([
        float(value.replace("D", "E").replace("d", "e"))
        for value in re.findall(r"[-+]?\d+(?:\.\d*)?(?:[DdEe][+-]?\d+)?", text)
    ])

def close(a, b, rtol=1.0e-9, atol=1.0e-10):
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b), 1.0)

fortran_params = parse_fortran_assignments(param_file)
fortran_init = parse_fortran_initial_values(init_file)
user_local_params = {"fS": 0.797, "fD": 1.0, "fCAo": 1.21}

param_mismatches = []
missing_from_fortran = []
for name, value in params.items():
    key = name.lower()
    source_value = user_local_params[name] if name in user_local_params else fortran_params.get(key)
    if source_value is None:
        missing_from_fortran.append(name)
    elif not close(value, source_value):
        param_mismatches.append((name, value, source_value))

initial_mismatches = [
    (f"x{i + 1}", float(x0[i]), float(fortran_init[i]))
    for i in range(31)
    if not close(float(x0[i]), float(fortran_init[i]))
]

def r0(name):
    return params[name] * P_CONV

def c0(name):
    return params[name] / P_CONV

def l0(name):
    return params[name] * P_CONV

expected_vessels = {
    "aorta_arch": {"R_poiseuille": 1.0, "C": params["fCAo"] * c0("Cao"), "L": 0.0},
    "upper_art": {"R_poiseuille": r0("Ruba"), "C": c0("Caub"), "L": l0("Luba")},
    "upper_cap": {"R_poiseuille": r0("Rub"), "C": c0("Cvub"), "L": 0.0},
    "upper_vein": {"R_poiseuille": r0("Rubv"), "C": 0.0, "L": 0.0},
    "svc": {"R_poiseuille": r0("Rsvc"), "C": c0("Csvc"), "L": l0("Luc")},
    "aorta_th": {"R_poiseuille": r0("Rthao"), "C": c0("Cthao"), "L": l0("Lthao")},
    "aorta_ab": {"R_poiseuille": r0("Rabao"), "C": c0("Cabao"), "L": l0("Labao")},
    "leg_art": {"R_poiseuille": r0("Rlega"), "C": c0("Clega"), "L": l0("Llega")},
    "leg_cap": {"R_poiseuille": r0("Rlegc"), "C": c0("Clegv"), "L": 0.0},
    "leg_vein": {"R_poiseuille": r0("Rlegv"), "C": c0("Cabivc"), "L": 0.0},
    "ivc_ab": {"R_poiseuille": r0("Rabivc"), "C": 0.0, "L": 0.0},
    "liver_art": {"R_poiseuille": r0("Rla"), "C": c0("Cl"), "L": 0.0},
    "liver_vein": {"R_poiseuille": r0("Rlv"), "C": 0.0, "L": 0.0},
    "kid_art": {"R_poiseuille": r0("Rka"), "C": c0("Ck"), "L": 0.0},
    "kid_vein": {"R_poiseuille": r0("Rkv"), "C": 0.0, "L": 0.0},
    "int_art": {"R_poiseuille": r0("Ria"), "C": c0("Ci"), "L": 0.0},
    "int_vein": {"R_poiseuille": r0("Riv"), "C": 0.0, "L": 0.0},
    "ivc": {"R_poiseuille": r0("Rthivc"), "C": c0("Cthivc"), "L": l0("Llc")},
    "Fontan": {"R_poiseuille": 1.0, "C": c0("Ctcpc"), "L": 0.0},
    "left_pul_art": {"R_poiseuille": r0("Rapoll"), "C": c0("Clungl"), "L": l0("Lapoll")},
    "left_pul_cap": {"R_poiseuille": r0("Rlungl"), "C": c0("Cvpl"), "L": l0("Llungl")},
    "left_pul_vein": {"R_poiseuille": r0("Rvpl"), "C": 0.0, "L": 0.0},
    "right_pul_art": {"R_poiseuille": r0("Rapolr"), "C": c0("Clungr"), "L": l0("Lapolr")},
    "right_pul_cap": {"R_poiseuille": r0("Rlungr"), "C": c0("Cvpr"), "L": l0("Llungr")},
    "right_pul_vein": {"R_poiseuille": r0("Rvpr"), "C": 0.0, "L": 0.0},
    "collateral_left": {"R_poiseuille": params["Rcoll"] / (1.0 - params["coll_split"]) * P_CONV, "C": 0.0, "L": 0.0},
    "collateral_right": {"R_poiseuille": params["Rcoll"] / params["coll_split"] * P_CONV, "C": 0.0, "L": 0.0},
}

vessel_mismatches = []
json_vessels = {v["vessel_name"]: v["zero_d_element_values"] for v in sv0d_config["vessels"]}
for vessel_name, expected_values in expected_vessels.items():
    actual_values = json_vessels[vessel_name]
    for key, expected_value in expected_values.items():
        if not close(actual_values[key], expected_value):
            vessel_mismatches.append((vessel_name, key, actual_values[key], expected_value))

def piecewise_activation_audit(t, activation):
    period = sv0d_config["simulation_parameters"].get("cardiac_period", 60.0 / params["HR"])
    phase = t % period
    phase = phase + period if phase < activation["contract_start"] else phase
    if phase < activation["contract_start"]:
        return 0.0
    if phase < activation["contract_start"] + activation["contract_duration"]:
        tau = (phase - activation["contract_start"]) / activation["contract_duration"]
        return 0.5 * (1.0 - np.cos(np.pi * tau))
    if phase < activation["relax_start"]:
        return 1.0
    if phase < activation["relax_start"] + activation["relax_duration"]:
        tau = (phase - activation["relax_start"]) / activation["relax_duration"]
        return 0.5 * (1.0 + np.cos(np.pi * tau))
    return 0.0

chambers_audit = {item["name"]: item for item in sv0d_config["chambers"]}
valves_audit = {item["name"]: item for item in sv0d_config["valves"]}

def sv0d_chamber_pressure(name, volume, t=0.0):
    chamber = chambers_audit[name]
    values = chamber["values"]
    act = piecewise_activation_audit(t, chamber["activation_function"])
    elastance = values["Epass"] + act * (values["Emax"] - values["Epass"])
    return elastance * (volume - values["Vrest"]) / P_CONV

def sv0d_valve_flow(upstream_pressure, downstream_pressure, name):
    valve = valves_audit[name]["params"]
    dp = (upstream_pressure - downstream_pressure) * P_CONV
    resistance = valve["Rmin"] if dp >= 0.0 else valve["Rmax"]
    return dp / resistance

p_la0 = sv0d_chamber_pressure("left_atrium", x0[26])
p_lv0 = sv0d_chamber_pressure("left_ventricle", x0[28])
q_mitral0 = sv0d_valve_flow(p_la0, p_lv0, "mitral")
q_aortic0 = sv0d_valve_flow(p_lv0, x0[0], "aortic")
q_la_in0 = (x0[21] - p_la0) / params["Rvpr"] + (x0[25] - p_la0) / params["Rvpl"]

expected_ic_subset = {
    "Vc:left_atrium": x0[26],
    "Vc:left_ventricle": x0[28],
    "flow:J2:left_atrium": q_la_in0,
    "flow:left_atrium:mitral": 0.0,
    "flow:mitral:left_ventricle": 0.0,
    "flow:left_ventricle:aortic": 0.0,
    "flow:aortic:aorta_arch": 0.0,
    "flow:ivc:J6": abs(x0[29]),
    "flow:svc:J6": abs(x0[30]),
    "pressure:J2:left_atrium": p_la0 * P_CONV,
    "pressure:left_atrium:mitral": p_la0 * P_CONV,
    "pressure:mitral:left_ventricle": p_lv0 * P_CONV,
    "pressure:left_ventricle:aortic": p_lv0 * P_CONV,
    "pressure:aortic:aorta_arch": x0[0] * P_CONV,
}
ic_mismatches = [
    (name, sv0d_ic.get(name), expected)
    for name, expected in expected_ic_subset.items()
    if name not in sv0d_ic or not close(float(sv0d_ic[name]), float(expected))
]

genbc_p_la0 = genbc["Psa"][0]
genbc_p_lv0 = genbc["Psv"][0]

def format_mismatches(rows, columns):
    if not rows:
        return "<p>none</p>"
    header = "".join(f"<th>{c}</th>" for c in columns)
    body = "".join("<tr>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>" for row in rows[:20])
    if len(rows) > 20:
        body += f"<tr><td colspan='{len(columns)}'>... {len(rows) - 20} more</td></tr>"
    return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"

summary_rows = [
    ("parameters_final.f scalar params", f"{len(params) - len(param_mismatches) - len(missing_from_fortran)}/{len(params)} exact", "fS/fD/fCAo are USER.f-local constants"),
    ("initial_values_final.f x1-x31", "match" if not initial_mismatches else f"{len(initial_mismatches)} mismatches", "x32-x37 are 3D-coupled face pressure states, not used in full 0D"),
    ("sv0D vessel R/C/L", "match" if not vessel_mismatches else f"{len(vessel_mismatches)} mismatches", "R,L multiplied by 1334; C divided by 1334"),
    ("sv0D valve/chamber port IC subset", "match" if not ic_mismatches else f"{len(ic_mismatches)} mismatches", "uses sv0D chamber pressure ICs and requested zero valve-flow ICs"),
]

intentional_rows = [
    ("chambers", "genBC nonlinear pressure-volume law is replaced by tuned LinearElastanceChamber"),
    ("valves", "genBC inertial/nonlinear valves are replaced by PiecewiseValve Rmin/Rmax"),
    ("Vrest", "genBC Vsa0/Vsv0 are not used as sv0D Vrest; current sv0D Vrest is 0.0"),
    ("IVC/SVC flow IC sign", "genBC x30/x31 start negative; sv0D stores positive flow in the vessel orientation into Fontan"),
    ("x32-x37", "Fortran 3D-coupled pulmonary face pressures are excluded from the standalone full 0D case"),
    ("aorta_arch and Fontan R", "dummy R_poiseuille=1.0 cgs is used because sv0D BloodVessel requires an R value"),
]

pressure_rows = [
    ("LA pressure at t=0", f"{genbc_p_la0:0.4g} mmHg", f"{p_la0:0.4g} mmHg", f"{p_la0 - genbc_p_la0:+0.4g}"),
    ("SV pressure at t=0", f"{genbc_p_lv0:0.4g} mmHg", f"{p_lv0:0.4g} mmHg", f"{p_lv0 - genbc_p_lv0:+0.4g}"),
    ("mitral flow implied at t=0", f"{x0[27]:0.4g} mL/s", f"{q_mitral0:0.4g} mL/s", "sv0D PiecewiseValve algebraic flow"),
    ("aortic flow implied at t=0", f"{genbc['Qao'][0]:0.4g} mL/s", f"{q_aortic0:0.4g} mL/s", "finite Rmax allows small backflow"),
]

display(HTML(
    "<h3>Parameter and initial-condition audit</h3>"
    + format_mismatches(summary_rows, ["check", "status", "note"])
    + "<h4>Current true mismatches</h4>"
    + format_mismatches(param_mismatches, ["name", "Python/sv0D", "Fortran"])
    + format_mismatches(initial_mismatches, ["state", "Python", "Fortran"])
    + format_mismatches(vessel_mismatches, ["vessel", "field", "JSON", "expected"])
    + format_mismatches(ic_mismatches, ["IC key", "JSON", "expected"])
    + "<h4>Intentional translations / not 1-to-1 with original genBC physics</h4>"
    + format_mismatches(intentional_rows, ["item", "reason"])
    + "<h4>Initial heart-law differences</h4>"
    + format_mismatches(pressure_rows, ["quantity", "original genBC", "sv0D law", "note"])
))
'''.strip()


PYTHON_SV0D_HEART_CODE = r'''
sv0d_config = json.loads(JSON_PATH.read_text())
chamber_config = {item["name"]: item for item in sv0d_config["chambers"]}
valve_config = {item["name"]: item for item in sv0d_config["valves"]}

def piecewise_activation(t, activation, p=params):
    period = sv0d_config["simulation_parameters"].get("cardiac_period", 60.0 / p["HR"])
    phase = t % period
    contract_start = activation["contract_start"]
    relax_start = activation["relax_start"]
    contract_duration = activation["contract_duration"]
    relax_duration = activation["relax_duration"]
    phase_unwrapped = phase + period if phase < contract_start else phase

    if phase_unwrapped < contract_start:
        return 0.0
    if phase_unwrapped < contract_start + contract_duration:
        tau = (phase_unwrapped - contract_start) / contract_duration
        return 0.5 * (1.0 - np.cos(np.pi * tau))
    if phase_unwrapped < relax_start:
        return 1.0
    if phase_unwrapped < relax_start + relax_duration:
        tau = (phase_unwrapped - relax_start) / relax_duration
        return 0.5 * (1.0 + np.cos(np.pi * tau))
    return 0.0

def linear_chamber_pressure_mmhg(t, volume, chamber):
    values = chamber["values"]
    act = piecewise_activation(t, chamber["activation_function"])
    elastance = values["Epass"] + act * (values["Emax"] - values["Epass"])
    return elastance * (volume - values["Vrest"]) / P_CONV

def piecewise_valve_flow_mls(p_up_mmhg, p_down_mmhg, valve):
    dp_cgs = (p_up_mmhg - p_down_mmhg) * P_CONV
    params_valve = valve["params"]
    resistance = params_valve["Rmin"] if dp_cgs >= 0.0 else params_valve["Rmax"]
    return dp_cgs / resistance

def rhs_python_sv0d_heart(x, p=params):
    f = np.zeros(32)
    y = {}
    t = x[31]
    tc = 60.0 / p["HR"]
    tr = 4.0 * tc

    psa = linear_chamber_pressure_mmhg(t, x[26], chamber_config["left_atrium"])
    psv = linear_chamber_pressure_mmhg(t, x[28], chamber_config["left_ventricle"])

    tmr = t % tr
    tsr = p["Trr"] * tr
    if tmr < tsr:
        dar = (p["pi"] / tsr) * np.sin(2.0 * p["pi"] * tmr / tsr)
    else:
        dar = 0.0
    dpith = p["Apith"] * dar

    qmitral = piecewise_valve_flow_mls(psa, psv, valve_config["mitral"])
    qao = piecewise_valve_flow_mls(psv, x[0], valve_config["aortic"])
    qub = (x[2] - x[3]) / p["Rub"]
    qvub = (x[3] - x[4]) / p["Rubv"]
    qsvc = x[30]
    qlegc = (x[10] - x[11]) / p["Rlegc"]
    qlegv = (abs(x[11] - x[12]) + x[11] - x[12]) / (2.0 * p["Rlegv"])
    qabivc = (x[12] - x[13]) / p["Rabivc"]
    qthivc = x[29]
    qla = (x[6] - x[14]) / p["Rla"]
    qlv = (x[14] - x[13]) / p["Rlv"]
    qia = (x[8] - x[15]) / p["Ria"]
    qiv = (x[15] - x[14]) / p["Riv"]
    qka = (x[6] - x[16]) / p["Rka"]
    qkv = (x[16] - x[13]) / p["Rkv"]
    qvpr = (x[21] - psa) / p["Rvpr"]
    qvpl = (x[25] - psa) / p["Rvpl"]
    qcoll_r = p["coll_split"] * (x[6] - x[21]) / p["Rcoll"]
    qcoll_l = (1.0 - p["coll_split"]) * (x[6] - x[25]) / p["Rcoll"]
    qcoll = qcoll_r + qcoll_l

    f[0] = (qao - x[1] - x[5]) / (p["fCAo"] * p["Cao"]) + dpith
    f[1] = (x[0] - x[1] * p["Ruba"] - x[2]) / p["Luba"]
    f[2] = (x[1] - qub) / p["Caub"]
    f[3] = (qub - qvub) / p["Cvub"]
    f[4] = (qvub - qsvc) / p["Csvc"] + dpith
    f[5] = (x[0] - x[5] * p["Rthao"] - x[6]) / p["Lthao"]
    f[6] = (x[5] - qcoll - x[7] - qla - qka) / p["Cthao"] + dpith
    f[7] = (x[6] - x[7] * p["Rabao"] - x[8]) / p["Labao"]
    f[8] = (x[7] - qia - x[9]) / p["Cabao"]
    f[9] = (x[8] - x[9] * p["Rlega"] - x[10]) / p["Llega"]
    f[10] = (x[9] - qlegc) / p["Clega"]
    f[11] = (qlegc - qlegv) / p["Clegv"]
    f[12] = (qlegv - qabivc) / p["Cabivc"]
    f[13] = (qabivc + qkv + qlv - qthivc) / p["Cthivc"] + dpith
    f[14] = (qiv + qla - qlv) / p["Cl"]
    f[15] = (qia - qiv) / p["Ci"]
    f[16] = (qka - qkv) / p["Ck"]
    f[17] = (qsvc + qthivc - x[18] - x[22]) / p["Ctcpc"] + dpith
    f[18] = (x[17] - x[18] * p["Rapolr"] - x[19]) / p["Lapolr"]
    f[19] = (x[18] - x[20]) / p["Clungr"] + dpith
    f[20] = (x[19] - x[20] * p["Rlungr"] - x[21]) / p["Llungr"]
    f[21] = (x[20] + qcoll_r - qvpr) / p["Cvpr"] + dpith
    f[22] = (x[17] - x[22] * p["Rapoll"] - x[23]) / p["Lapoll"]
    f[23] = (x[22] - x[24]) / p["Clungl"] + dpith
    f[24] = (x[23] - x[24] * p["Rlungl"] - x[25]) / p["Llungl"]
    f[25] = (x[24] + qcoll_l - qvpl) / p["Cvpl"] + dpith
    f[26] = qvpr + qvpl - qmitral
    f[27] = 0.0
    f[28] = qmitral - qao
    f[29] = (x[13] - x[29] * p["Rthivc"] - x[17]) / p["Llc"]
    f[30] = (x[4] - x[30] * p["Rsvc"] - x[17]) / p["Luc"]
    f[31] = 1.0

    y.update(
        Psa=psa,
        Psv=psv,
        Qao=qao,
        Qmitral=qmitral,
        Qcoll_l=qcoll_l,
        Qcoll_r=qcoll_r,
        Qcoll=qcoll,
    )
    return f, y

def run_python_sv0d_heart_reference():
    steps = 250
    substeps = 40
    ncyc = 20
    dt = (60.0 / params["HR"]) / (steps * substeps)
    total_steps = ncyc * steps * substeps
    out = []
    x = x0.copy()

    for i in range(total_steps + 1):
        if i % substeps == 0:
            _, y = rhs_python_sv0d_heart(x)
            if i == 0:
                y = {**y, "Qao": 0.0, "Qmitral": 0.0}
            out.append([
                x[31], *x[:31], y["Psa"], y["Psv"], y["Qao"],
                y["Qmitral"], y["Qcoll_l"], y["Qcoll_r"], y["Qcoll"]
            ])
        if i == total_steps:
            break
        k1, _ = rhs_python_sv0d_heart(x)
        k2, _ = rhs_python_sv0d_heart(x + 0.5 * dt * k1)
        k3, _ = rhs_python_sv0d_heart(x + 0.5 * dt * k2)
        k4, _ = rhs_python_sv0d_heart(x + dt * k3)
        x = x + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

    cols = ["time"] + [f"x{i}" for i in range(1, 32)] + ["Psa", "Psv", "Qao", "Qmitral", "Qcoll_l", "Qcoll_r", "Qcoll"]
    out_path = CASE_DIR / "python_genbc_sv0d_heart_result.csv"
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(cols)
        writer.writerows(out)
    arr = np.asarray(out, dtype=float)
    return {name: arr[:, i] for i, name in enumerate(cols)}, out_path

hybrid, hybrid_out = run_python_sv0d_heart_reference()

def hybrid_col(name):
    return hybrid[name]

hybrid_compare = {
    "time": compare["time"],
    "hybrid_Pao": hybrid_col("x1"),
    "hybrid_Pfontan": hybrid_col("x18"),
    "hybrid_Psa": hybrid_col("Psa"),
    "hybrid_Psv": hybrid_col("Psv"),
    "hybrid_Qao": hybrid_col("Qao"),
    "hybrid_Qmitral": hybrid_col("Qmitral"),
    "hybrid_Vla": hybrid_col("x27"),
    "hybrid_Vlv": hybrid_col("x29"),
    "hybrid_Qcollateral": hybrid_col("Qcoll"),
    "hybrid_Qlpa": hybrid_col("x23"),
    "hybrid_Qrpa": hybrid_col("x19"),
    "hybrid_Qivc": hybrid_col("x30"),
    "hybrid_Qsvc": hybrid_col("x31"),
}

last = compare["time"] >= (compare["time"].max() - 60.0 / params["HR"])
all_compare = {**compare, **hybrid_compare}

def rel_l2_arrays(a, b):
    aa = np.asarray(a, dtype=float)[last]
    bb = np.asarray(b, dtype=float)[last]
    return float(np.linalg.norm(aa - bb) / (np.linalg.norm(bb) + 1.0e-12))

def metric_value(name, reducer=np.mean):
    return float(reducer(all_compare[name][last]))

def beat_integral(name):
    return float(np.trapezoid(all_compare[name][last], all_compare["time"][last]))

def volume_swing(name):
    values = all_compare[name][last]
    return float(np.max(values) - np.min(values))

def metric_row(label, genbc_name, hybrid_name, sv0d_name, reducer=np.mean):
    return [
        label,
        metric_value(genbc_name, reducer),
        metric_value(hybrid_name, reducer),
        metric_value(sv0d_name, reducer),
        rel_l2_arrays(all_compare[sv0d_name], all_compare[genbc_name]),
        rel_l2_arrays(all_compare[sv0d_name], all_compare[hybrid_name]),
    ]

def integral_row(label, genbc_name, hybrid_name, sv0d_name):
    return [
        label,
        beat_integral(genbc_name),
        beat_integral(hybrid_name),
        beat_integral(sv0d_name),
        rel_l2_arrays(all_compare[sv0d_name], all_compare[genbc_name]),
        rel_l2_arrays(all_compare[sv0d_name], all_compare[hybrid_name]),
    ]

def swing_row(label, genbc_name, hybrid_name, sv0d_name):
    return [
        label,
        volume_swing(genbc_name),
        volume_swing(hybrid_name),
        volume_swing(sv0d_name),
        rel_l2_arrays(all_compare[sv0d_name], all_compare[genbc_name]),
        rel_l2_arrays(all_compare[sv0d_name], all_compare[hybrid_name]),
    ]

three_way_slide_metrics = [
    ["Pressure", "Pao", "Mean", "mmHg", metric_value("genBC_Pao"), metric_value("hybrid_Pao"), metric_value("sv0D_Pao")],
    ["Pressure", "Pao", "Max", "mmHg", metric_value("genBC_Pao", np.max), metric_value("hybrid_Pao", np.max), metric_value("sv0D_Pao", np.max)],
    ["Pressure", "Fontan", "Mean", "mmHg", metric_value("genBC_Pfontan"), metric_value("hybrid_Pfontan"), metric_value("sv0D_Pfontan")],
    ["Pressure", "Psa / Pla", "Mean", "mmHg", metric_value("genBC_Pla"), metric_value("hybrid_Psa"), metric_value("sv0D_Pla")],
    ["Pressure", "Psa / Pla", "Max", "mmHg", metric_value("genBC_Pla", np.max), metric_value("hybrid_Psa", np.max), metric_value("sv0D_Pla", np.max)],
    ["Pressure", "Psv / Plv", "Mean", "mmHg", metric_value("genBC_Plv"), metric_value("hybrid_Psv"), metric_value("sv0D_Plv")],
    ["Pressure", "Psv / Plv", "Max", "mmHg", metric_value("genBC_Plv", np.max), metric_value("hybrid_Psv", np.max), metric_value("sv0D_Plv", np.max)],
    ["Flow", "Qao", "Mean", "mL/s", metric_value("genBC_Qao"), metric_value("hybrid_Qao"), metric_value("sv0D_Qao")],
    ["Flow", "Qao", "Max", "mL/s", metric_value("genBC_Qao", np.max), metric_value("hybrid_Qao", np.max), metric_value("sv0D_Qao", np.max)],
    ["Flow", "Collateral", "Mean", "mL/s", metric_value("genBC_Qcollateral"), metric_value("hybrid_Qcollateral"), metric_value("sv0D_Qcollateral")],
    ["Flow", "Collateral", "Max", "mL/s", metric_value("genBC_Qcollateral", np.max), metric_value("hybrid_Qcollateral", np.max), metric_value("sv0D_Qcollateral", np.max)],
    ["Flow", "LPA", "Mean", "mL/s", metric_value("genBC_Qlpa"), metric_value("hybrid_Qlpa"), metric_value("sv0D_Qlpa")],
    ["Flow", "RPA", "Mean", "mL/s", metric_value("genBC_Qrpa"), metric_value("hybrid_Qrpa"), metric_value("sv0D_Qrpa")],
    ["Flow", "IVC", "Mean", "mL/s", metric_value("genBC_Qivc"), metric_value("hybrid_Qivc"), metric_value("sv0D_Qivc")],
    ["Flow", "SVC", "Mean", "mL/s", metric_value("genBC_Qsvc"), metric_value("hybrid_Qsvc"), metric_value("sv0D_Qsvc")],
]

def rel_error(value, baseline):
    return abs(value - baseline) / (abs(baseline) + 1.0e-12)

def three_way_fmt_value(value):
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 10:
        return f"{value:.2f}"
    return f"{value:.3f}"

def three_way_fmt_pct(value):
    return f"{100.0 * value:.1f}%"

def three_way_error_class(err):
    if err < 0.10:
        return "good"
    if err < 0.25:
        return "moderate"
    return "large"

three_way_rows = []
for group, variable, statistic, unit, genbc_value, hybrid_value, sv0d_value in three_way_slide_metrics:
    three_way_rows.append({
        "group": group,
        "metric": f"{variable} {statistic.lower()}",
        "unit": unit,
        "genBC": genbc_value,
        "Python hybrid": hybrid_value,
        "sv0D": sv0d_value,
        "hybrid_rel_error": rel_error(hybrid_value, genbc_value),
        "sv0d_rel_error": rel_error(sv0d_value, genbc_value),
    })

three_way_csv = CASE_DIR / "fontan_0d_three_way_summary_slide_table.csv"
with three_way_csv.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["group", "metric", "unit", "genBC", "Python hybrid", "sv0D", "hybrid rel error", "sv0D rel error"])
    for row in three_way_rows:
        writer.writerow([
            row["group"], row["metric"], row["unit"], row["genBC"], row["Python hybrid"], row["sv0D"],
            row["hybrid_rel_error"], row["sv0d_rel_error"],
        ])

three_way_html_rows = "\n".join(
    "<tr>"
    f"<td class='group {row['group'].lower()}'>{row['group']}</td>"
    f"<td class='metric'>{row['metric']}</td>"
    f"<td class='unit'>{row['unit']}</td>"
    f"<td>{three_way_fmt_value(row['genBC'])}</td>"
    f"<td>{three_way_fmt_value(row['Python hybrid'])}</td>"
    f"<td>{three_way_fmt_value(row['sv0D'])}</td>"
    f"<td><span class='badge {three_way_error_class(row['hybrid_rel_error'])}'>{three_way_fmt_pct(row['hybrid_rel_error'])}</span></td>"
    f"<td><span class='badge {three_way_error_class(row['sv0d_rel_error'])}'>{three_way_fmt_pct(row['sv0d_rel_error'])}</span></td>"
    "</tr>"
    for row in three_way_rows
)

three_way_table_html = f"""
<style>
.three-way-slide-card {{
  width: 1040px;
  background: #ffffff;
  border: 1px solid #d8dee8;
  border-radius: 8px;
  overflow: hidden;
  font-family: Arial, Helvetica, sans-serif;
  color: #172033;
}}
.three-way-slide-title {{
  padding: 14px 18px 8px 18px;
  font-size: 20px;
  font-weight: 700;
}}
.three-way-slide-subtitle {{
  padding: 0 18px 12px 18px;
  color: #5b6576;
  font-size: 12px;
}}
.three-way-slide-card table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}}
.three-way-slide-card th {{
  background: #eef2f7;
  color: #334155;
  padding: 8px 10px;
  border-top: 1px solid #d8dee8;
  border-bottom: 1px solid #d8dee8;
  text-align: right;
  font-weight: 700;
}}
.three-way-slide-card th:nth-child(1),
.three-way-slide-card th:nth-child(2),
.three-way-slide-card th:nth-child(3) {{
  text-align: left;
}}
.three-way-slide-card td {{
  padding: 7px 10px;
  border-bottom: 1px solid #edf1f6;
  text-align: right;
  white-space: nowrap;
}}
.three-way-slide-card tr:nth-child(even) td {{
  background: #fafbfd;
}}
.three-way-slide-card td.metric,
.three-way-slide-card td.unit,
.three-way-slide-card td.group {{
  text-align: left;
}}
.three-way-slide-card td.unit {{
  color: #64748b;
}}
.three-way-slide-card td.group {{
  font-weight: 700;
  border-left: 5px solid #2b6cb0;
}}
.three-way-slide-card td.group.flow {{
  border-left-color: #0f8b6c;
}}
.three-way-slide-card .badge {{
  display: inline-block;
  min-width: 52px;
  padding: 3px 7px;
  border-radius: 999px;
  text-align: center;
  font-weight: 700;
}}
.three-way-slide-card .badge.good {{
  background: #e7f6ef;
  color: #087049;
}}
.three-way-slide-card .badge.moderate {{
  background: #fff3d8;
  color: #9a5a00;
}}
.three-way-slide-card .badge.large {{
  background: #fde8e8;
  color: #b42318;
}}
</style>
<div class="three-way-slide-card">
  <div class="three-way-slide-title">Full 0D Last-Cycle Comparison</div>
  <div class="three-way-slide-subtitle">genBC baseline, Python hybrid, and sv0D. Error columns are scalar relative errors versus genBC for each displayed mean or max value.</div>
  <table>
    <thead>
      <tr><th>Group</th><th>Metric</th><th>Unit</th><th>genBC</th><th>Python hybrid</th><th>sv0D</th><th>Hybrid error</th><th>sv0D error</th></tr>
    </thead>
    <tbody>{three_way_html_rows}</tbody>
  </table>
</div>
"""

display(HTML(
    "<h3>Three-way comparison: genBC, Python hybrid, and sv0D</h3>"
    f"<p>Wrote <code>{hybrid_out}</code>. The Python hybrid uses the genBC vascular ODEs with the sv0D chamber and valve laws; sv0D uses the full block/junction DAE formulation from the JSON.</p>"
    + three_way_table_html
))

svg_width = 1320
row_h = 34
top = 94
left = 28
columns = [left, 130, 380, 490, 625, 790, 925, 1090, 1260]
svg_height = top + row_h * (len(three_way_rows) + 1) + 34
parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_width}" height="{svg_height}" viewBox="0 0 {svg_width} {svg_height}">',
    '<rect width="100%" height="100%" fill="#ffffff"/>',
    '<rect x="14" y="14" width="1292" height="{0}" rx="8" fill="#ffffff" stroke="#d8dee8"/>'.format(svg_height - 28),
    svg_text(34, 48, "Full 0D Last-Cycle Comparison", 24, "700"),
    svg_text(34, 72, "genBC baseline, Python hybrid, and sv0D; errors = scalar relative errors versus genBC", 13, "400", fill="#5b6576"),
    f'<rect x="{left}" y="{top}" width="{columns[-1] - left}" height="{row_h}" fill="#eef2f7"/>',
]
headers = ["Group", "Metric", "Unit", "genBC", "Python hybrid", "sv0D", "Hybrid error", "sv0D error"]
for i, header in enumerate(headers):
    anchor = "start" if i < 3 else "end"
    x = columns[i] + 10 if i < 3 else columns[i + 1] - 10
    parts.append(svg_text(x, top + 23, header, 13, "700", anchor, "#334155"))

for idx, row in enumerate(three_way_rows):
    y = top + row_h * (idx + 1)
    fill = "#fafbfd" if idx % 2 else "#ffffff"
    accent = "#2b6cb0" if row["group"] == "Pressure" else "#0f8b6c"
    hybrid_badge_fill, hybrid_badge_text = {
        "good": ("#e7f6ef", "#087049"),
        "moderate": ("#fff3d8", "#9a5a00"),
        "large": ("#fde8e8", "#b42318"),
    }[three_way_error_class(row["hybrid_rel_error"])]
    sv0d_badge_fill, sv0d_badge_text = {
        "good": ("#e7f6ef", "#087049"),
        "moderate": ("#fff3d8", "#9a5a00"),
        "large": ("#fde8e8", "#b42318"),
    }[three_way_error_class(row["sv0d_rel_error"])]
    parts.extend([
        f'<rect x="{left}" y="{y}" width="{columns[-1] - left}" height="{row_h}" fill="{fill}"/>',
        f'<rect x="{left}" y="{y}" width="5" height="{row_h}" fill="{accent}"/>',
        svg_text(columns[0] + 10, y + 22, row["group"], 13, "700"),
        svg_text(columns[1] + 10, y + 22, row["metric"], 13),
        svg_text(columns[2] + 10, y + 22, row["unit"], 13, fill="#64748b"),
        svg_text(columns[4] - 10, y + 22, three_way_fmt_value(row["genBC"]), 13, anchor="end"),
        svg_text(columns[5] - 10, y + 22, three_way_fmt_value(row["Python hybrid"]), 13, anchor="end"),
        svg_text(columns[6] - 10, y + 22, three_way_fmt_value(row["sv0D"]), 13, anchor="end"),
        f'<rect x="{columns[6] + 20}" y="{y + 7}" width="76" height="20" rx="10" fill="{hybrid_badge_fill}"/>',
        svg_text(columns[7] - 20, y + 22, three_way_fmt_pct(row["hybrid_rel_error"]), 13, "700", "end", hybrid_badge_text),
        f'<rect x="{columns[7] + 20}" y="{y + 7}" width="76" height="20" rx="10" fill="{sv0d_badge_fill}"/>',
        svg_text(columns[8] - 20, y + 22, three_way_fmt_pct(row["sv0d_rel_error"]), 13, "700", "end", sv0d_badge_text),
    ])

parts.append("</svg>")
three_way_svg_path = CASE_DIR / "fontan_0d_three_way_summary_slide_table.svg"
three_way_svg_path.write_text("\n".join(parts))
display(HTML(f"<p>Wrote slide-ready three-way table: <code>{three_way_svg_path}</code> and <code>{three_way_csv}</code></p>"))

plot_window = compare["time"] >= (compare["time"].max() - 2 * (60.0 / params["HR"]))
pv_window = compare["time"] >= (compare["time"].max() - 60.0 / params["HR"])

def line_points(x, y, xmin, xmax, ymin, ymax, width, height, pad):
    xs = pad + (x - xmin) / max(xmax - xmin, 1.0e-12) * (width - 2 * pad)
    ys = height - pad - (y - ymin) / max(ymax - ymin, 1.0e-12) * (height - 2 * pad)
    return " ".join(f"{xx:0.1f},{yy:0.1f}" for xx, yy in zip(xs, ys))

def svg_three_way_timeseries(title, ylabel, genbc_name, hybrid_name, sv0d_name, width=420, height=250):
    pad = 42
    x = all_compare["time"][plot_window]
    series = [
        ("genBC", all_compare[genbc_name][plot_window], "#1f77b4", ""),
        ("Python hybrid", all_compare[hybrid_name][plot_window], "#2ca02c", ""),
        ("sv0D", all_compare[sv0d_name][plot_window], "#d62728", "6 4"),
    ]
    ymin = float(min(np.min(y) for _, y, _, _ in series))
    ymax = float(max(np.max(y) for _, y, _, _ in series))
    margin = 0.08 * max(ymax - ymin, 1.0)
    ymin -= margin
    ymax += margin
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    polylines = []
    for label, y, color, dash in series:
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        polylines.append(
            f'<polyline points="{line_points(x, y, xmin, xmax, ymin, ymax, width, height, pad)}" fill="none" stroke="{color}" stroke-width="2"{dash_attr}/>'
        )
    legend = (
        f'<text x="{pad+8}" y="{height-pad-12}" font-size="11" font-family="sans-serif" fill="#1f77b4">genBC</text>'
        f'<text x="{pad+68}" y="{height-pad-12}" font-size="11" font-family="sans-serif" fill="#2ca02c">Python hybrid</text>'
        f'<text x="{pad+186}" y="{height-pad-12}" font-size="11" font-family="sans-serif" fill="#d62728">sv0D</text>'
    )
    return f"""
    <svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" style="max-width:100%;border:1px solid #ddd;background:white">
      <text x="{pad}" y="22" font-size="15" font-family="sans-serif">{title}</text>
      <line x1="{pad}" y1="{height-pad}" x2="{width-pad}" y2="{height-pad}" stroke="#555"/>
      <line x1="{pad}" y1="{pad}" x2="{pad}" y2="{height-pad}" stroke="#555"/>
      <text x="{width-pad-80}" y="{height-12}" font-size="11" font-family="sans-serif">time [s]</text>
      <text x="8" y="{pad-10}" font-size="11" font-family="sans-serif">{ylabel}</text>
      {''.join(polylines)}
      {legend}
    </svg>
    """

def svg_three_way_pv_loop(title, genbc_volume, genbc_pressure, hybrid_volume, hybrid_pressure, sv0d_volume, sv0d_pressure, width=420, height=250):
    pad = 42
    series = [
        ("genBC", all_compare[genbc_volume][pv_window], all_compare[genbc_pressure][pv_window], "#1f77b4", ""),
        ("Python hybrid", all_compare[hybrid_volume][pv_window], all_compare[hybrid_pressure][pv_window], "#2ca02c", ""),
        ("sv0D", all_compare[sv0d_volume][pv_window], all_compare[sv0d_pressure][pv_window], "#d62728", "6 4"),
    ]
    xmin = float(min(np.min(x) for _, x, _, _, _ in series))
    xmax = float(max(np.max(x) for _, x, _, _, _ in series))
    ymin = float(min(np.min(y) for _, _, y, _, _ in series))
    ymax = float(max(np.max(y) for _, _, y, _, _ in series))
    xmargin = 0.08 * max(xmax - xmin, 1.0)
    ymargin = 0.08 * max(ymax - ymin, 1.0)
    xmin -= xmargin
    xmax += xmargin
    ymin -= ymargin
    ymax += ymargin
    polylines = []
    for label, x, y, color, dash in series:
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        polylines.append(
            f'<polyline points="{line_points(x, y, xmin, xmax, ymin, ymax, width, height, pad)}" fill="none" stroke="{color}" stroke-width="2"{dash_attr}/>'
        )
    legend = (
        f'<text x="{pad+8}" y="{height-pad-12}" font-size="11" font-family="sans-serif" fill="#1f77b4">genBC</text>'
        f'<text x="{pad+68}" y="{height-pad-12}" font-size="11" font-family="sans-serif" fill="#2ca02c">Python hybrid</text>'
        f'<text x="{pad+186}" y="{height-pad-12}" font-size="11" font-family="sans-serif" fill="#d62728">sv0D</text>'
    )
    return f"""
    <svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" style="max-width:100%;border:1px solid #ddd;background:white">
      <text x="{pad}" y="22" font-size="15" font-family="sans-serif">{title}</text>
      <line x1="{pad}" y1="{height-pad}" x2="{width-pad}" y2="{height-pad}" stroke="#555"/>
      <line x1="{pad}" y1="{pad}" x2="{pad}" y2="{height-pad}" stroke="#555"/>
      <text x="{width-pad-92}" y="{height-12}" font-size="11" font-family="sans-serif">Volume [mL]</text>
      <text x="8" y="{pad-10}" font-size="11" font-family="sans-serif">Pressure [mmHg]</text>
      {''.join(polylines)}
      {legend}
    </svg>
    """

three_way_plots = [
    svg_three_way_timeseries("Pao", "Pressure [mmHg]", "genBC_Pao", "hybrid_Pao", "sv0D_Pao"),
    svg_three_way_timeseries("Psv / Plv", "Pressure [mmHg]", "genBC_Plv", "hybrid_Psv", "sv0D_Plv"),
    svg_three_way_timeseries("Qao", "Flow [mL/s]", "genBC_Qao", "hybrid_Qao", "sv0D_Qao"),
    svg_three_way_timeseries("Mitral flow", "Flow [mL/s]", "genBC_Qmitral", "hybrid_Qmitral", "sv0D_Qmitral"),
]
display(HTML("<div style='display:grid;grid-template-columns:repeat(2,minmax(320px,1fr));gap:14px'>" + "".join(three_way_plots) + "</div>"))

three_way_pv_plots = [
    svg_three_way_pv_loop("Left atrium PV loop", "genBC_Vla", "genBC_Pla", "hybrid_Vla", "hybrid_Psa", "sv0D_Vla", "sv0D_Pla"),
    svg_three_way_pv_loop("Single ventricle PV loop", "genBC_Vlv", "genBC_Plv", "hybrid_Vlv", "hybrid_Psv", "sv0D_Vlv", "sv0D_Plv"),
]
display(HTML("<div style='display:grid;grid-template-columns:repeat(2,minmax(320px,1fr));gap:14px;margin-top:14px'>" + "".join(three_way_pv_plots) + "</div>"))
'''.strip()


def build_notebook() -> dict:
    return {
        "cells": [
            markdown_cell(
                "# genBC vs sv0D Fontan LPN comparison\n\n"
                "This notebook runs the closed-loop Python version of the old genBC LPN, runs the translated svZeroDSolver JSON, and compares last-cycle waveforms/metrics. "
                "The vascular topology and vascular parameters follow the genBC model directly. The chamber/valve layer uses the older explicit-heart svZeroDSolver schema: `LinearElastanceChamber` and `PiecewiseValve`."
            ),
            code_cell(NOTEBOOK_CODE),
            code_cell(GENBC_MODEL_CODE),
            code_cell(RUN_AND_COMPARE_CODE),
            code_cell(METRICS_AND_PLOTS_CODE),
            code_cell(VALVE_ANALYSIS_CODE),
            code_cell(PARAMETER_IC_AUDIT_CODE),
            code_cell(PYTHON_SV0D_HEART_CODE),
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(build_config(), indent=2) + "\n")
    NOTEBOOK_PATH.write_text(json.dumps(build_notebook(), indent=2) + "\n")
    print(f"Wrote {JSON_PATH}")
    print(f"Wrote {NOTEBOOK_PATH}")


if __name__ == "__main__":
    main()
