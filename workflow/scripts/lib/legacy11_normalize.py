#!/usr/bin/env python3
# =============================================================================
# legacy11_normalize.py  -  Stage 7 label normalization (vendored)
# =============================================================================
# VENDORED VERBATIM from the legacy Stage 7 script:
#     context/Legacy_stage7/11_final_reordenar_y_normalizar_v6.py
#
# Everything below the header is an unmodified copy of that script's
# label-normalization layer (from the top of the file down to, but excluding,
# `def dataset_from_name`, which belongs to the legacy I/O layer that
# Stage 7 replaces).
#
# WHY VENDORED RATHER THAN REWRITTEN. This layer encodes the published
# ConotoxinBenchmark methodology (curated exact mappings, rule-based
# transformations and the confidence decision tree). Any paraphrase risks
# silently changing published results, so the code is copied byte-for-byte
# and only its I/O wrapper is reimplemented (scripts 31 / 32).
#
# The only edit applied is removal of `from __future__ import annotations`
# (unnecessary on the pinned python=3.11 and a hazard if this module is ever
# executed through Snakemake's `script:` directive -- see gotcha #1).
#
# DO NOT EDIT the logic here. To change normalization or classification
# behaviour, open an Issue and amend the legacy source + this vendored copy
# together, so provenance stays auditable.
# =============================================================================

# -*- coding: utf-8 -*-

import argparse
import math
import re
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

FINAL_ORDER = [
    "encabezados",
    "conoserver",
    "ncbi",
    "uniprot",
    "rnabloom_spades_trinity_regex_Sequence Frequency",
    "rnabloom_spades_trinity_phmm_Sequence Frequency",
    "rnabloom_spades_trinity_phmm_E-value Superfamily",
    "# A.A",
    "# Cysteine(s)",
    "% Hydrophobicity (Signal)",
    "Superfamily (Signal)",
    "Protein Sequence",
]

MOVE_TO_MOVIDAS = [
    "Superfamily (Mature)",
    "Superfamily (Pro-region)",
]

WIDTHS_PX = [165, 102, 56, 202, 265, 265, 350, 31, 65, 125, 104, 1242]
WIDTHS_XL = [round((px - 5) / 7, 2) for px in WIDTHS_PX]

EMPTY_TOKENS = {"", "nan", "none", "null", "na", "n/a", "-", "—"}

TO_EMPTY = {"PVIC", "PVIB", "PVIE"}

FINAL_SUPERFAMILIES = {
    "A","B","B1","B2","B3","C","D","E","E1","F","G","G2","H","I","I1","I2","I3","J","K","L","M","N",
    "O","O1","O2","O3","P","Q","R","S","T","U","V","W","Y","Z"
}

FINAL_SPECIAL = {
    "Conantokin","Conkunitzin","Conopressin","Contulakin",
    "con-ikot-ikot","Conodipine","Conomarphin","Contryphan",
    "Conoporin","Conorfamide","Insulin",
    "Consomatin","ConoCAP","Elevenin","Scratcher",
    "NO_conotoxin","Conotoxina_no_clasificada"
}

# Familias divergentes válidas de ConoSorter/ConoServer.
# Se elimina únicamente el prefijo "Divergent"; el motivo se conserva.
DIVERGENT_SUPERFAMILIES = {
    "MSTLGMTLL-",
    "M---L-LTVA",
    "MKFPLLFISL",
    "MKLCVVIVLL",
    "MSKLVILAVL",
    "MRFLHFLIVA",
    "MRFYIGLMAA",
}

# Etiquetas válidas observadas directamente en Superfamily (Signal).
# Se conservan como categorías independientes; no se homologan mediante
# ConoServer, NCBI o UniProt.
SIGNAL_ONLY_FAMILIES = {
    "M2", "O4", "Y2", "Y3", "H2", "N2", "X",
    "MTFLLLLVSV", "MRCLSIFVLL", "MKLLLTLLLG",
    "MKVAVVLLVS", "MTAKATLLVL", "MTLTFLLVVA",
}

VALID_FINAL = (
    FINAL_SUPERFAMILIES
    | FINAL_SPECIAL
    | DIVERGENT_SUPERFAMILIES
    | SIGNAL_ONLY_FAMILIES
)

def normalize_space(x: object) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and math.isnan(x):
        return ""
    s = str(x)
    # remove zero-width chars
    s = re.sub(r'[\u200B-\u200D\uFEFF]', '', s)
    # normalize unicode dashes
    s = re.sub(r'[\u2010\u2011\u2012\u2013\u2014\u2015\u2212]', '-', s)
    s = re.sub(r'\s+', ' ', s)
    return s.strip()

def canon_key(x: object) -> str:
    s = normalize_space(x)
    # normalize spaces around punctuation that tend to vary in Excel
    s = re.sub(r'\s+\)', ')', s)
    s = re.sub(r'\(\s+', '(', s)
    s = re.sub(r'\s+,', ',', s)
    s = re.sub(r',\s+', ', ', s)
    s = re.sub(r'\s+_', '_', s)
    s = re.sub(r'_\s+', '_', s)
    s = re.sub(r'\s*-\s*', '-', s)
    s = re.sub(r'\s+', ' ', s).strip()
    return s.casefold()

def is_empty_token(s: str) -> bool:
    return canon_key(s) in {t.casefold() for t in EMPTY_TOKENS}

def clean_signal_value(s: object) -> str:
    t = normalize_space(s)
    if not t:
        return ""
    t = re.sub(r'^Divergent\s+', '', t, flags=re.I).strip()
    t = re.sub(r'\s*\((?:Signal|Pro-region)\)\s*$', '', t, flags=re.I).strip()
    if re.match(r'^SF-mi\d+$', t, flags=re.I):
        return "Conotoxina_no_clasificada"
    return t

FAMILY_ALT = r'(?:B1|B2|B3|I1|I2|I3|O1|O2|O3|G2|A|B|C|D|E|F|G|H|I|J|K|L|M|N|O|P|Q|R|S|T|U|V|W|Y|Z)'

def clean_conoserver_value(s: object) -> str:
    """
    Normaliza exclusivamente la primera anotación ConoServer.

    El script 8 concatena los valores en el orden del TSV DIAMOND. DIAMOND,
    sin --top, reporta los targets ordenados por e-value; por tanto se conserva
    el primer valor histórico y no se usa ninguna otra fuente para decidirlo.
    """
    t = normalize_space(s)
    if not t:
        return ""

    # Conservar solo la anotación del primer hit reportado.
    first = normalize_space(t.split(",", 1)[0])
    if not first:
        return ""

    # Las familias Divergent son etiquetas válidas e independientes.
    first = re.sub(r'^Divergent\s+', '', first, flags=re.I).strip()

    # Homologar "<etiqueta> superfamily" a "<etiqueta>".
    first = re.sub(r'\s+superfamily\s*$', '', first, flags=re.I).strip()

    special = {
        "insulin": "Insulin",
        "conkunitzin": "Conkunitzin",
        "conantokin": "Conantokin",
        "conopressin": "Conopressin",
        "contulakin": "Contulakin",
        "con-ikot-ikot": "con-ikot-ikot",
        "conodipine": "Conodipine",
        "conomarphin": "Conomarphin",
        "contryphan": "Contryphan",
        "conoporin": "Conoporin",
        "conorfamide": "Conorfamide",
    }
    key = canon_key(first)
    if key in special:
        return special[key]

    # Familias alfanuméricas canónicas.
    upper = first.upper()
    if upper in FINAL_SUPERFAMILIES:
        return upper

    # Motivos Divergent: conservar mayúsculas y guiones tal como están.
    if first in DIVERGENT_SUPERFAMILIES:
        return first

    return first


EXACT_MAP = {
    # direct family/class mappings
    "A2": "A",
    "B-": "B1",
    "Linear-P": "B2",
    "Conantokin E1": "E1",
    "J-superfamily Vt14.2": "J",
    "M-superfamily": "M",
    "G066 VD Superfamily B": "B",
    "G065 VD Superfamily B": "B",
    "G129 RS Superfamily M": "M",
    "G130 RS Superfamily O1": "O1",
    "Bromosleeper O3-superfamily protein": "O3",

    # peptides / classes
    "Conophysin-conopressin": "Conopressin",
    "Conopressin/conophysin": "Conopressin",
    "Hormone superfamily conopressin/conophysin": "Conopressin",
    "contulakin": "Contulakin",
    "Con-ikot-ikot_Vc21.1 prepropeptide": "con-ikot-ikot",

    # conantokin / conkunitzin / contulakin exacts
    "G006_VD_Conantokin__": "Conantokin",
    "G021_VD_Conantokin__": "Conantokin",
    "G026_VD_Conantokin__": "Conantokin",
    "G027_VD_Conantokin__": "Conantokin",
    "G035_VD_Conantokin__": "Conantokin",
    "G037_VD_Conkunitzin__": "Conkunitzin",
    "G038_VD_Conkunitzin__": "Conkunitzin",
    "G041_VD_Contulakin__": "Contulakin",

    # Conodipine / Conomarphin / Contryphan / Conorfamide family exacts
    "Conodipine 5": "Conodipine",
    "Conodipine_Vc2 prepropeptide": "Conodipine",
    "Conodipine_Vc3 prepropeptide": "Conodipine",
    "Conomarphin analog Bt1": "Conomarphin",
    "Conomarphin analog Bt2": "Conomarphin",
    "Conomarphin analog Mr1": "Conomarphin",
    "Conomarphin_Mr1__1": "Conomarphin",
    "Conomarphin_Mr1__3": "Conomarphin",
    "Conomarphin_Mr2_": "Conomarphin",
    "Contryphan_M__1": "Contryphan",

    # unclassified families / tags
    "Cerm-01": "Conotoxina_no_clasificada",
    "Cerm-06": "Conotoxina_no_clasificada",
    "Cerm-07": "Conotoxina_no_clasificada",
    "Cerm-10": "Conotoxina_no_clasificada",
    "Cerm-11": "Conotoxina_no_clasificada",
    "Cerm-13": "Conotoxina_no_clasificada",
    "Ggeo03-1": "Conotoxina_no_clasificada",
    "Ggeo03-2": "Conotoxina_no_clasificada",
    "Rimp-05": "Conotoxina_no_clasificada",
    "Superfamily Rimp-01": "Conotoxina_no_clasificada",
    "Superfamily Pmag-02": "Conotoxina_no_clasificada",
    "Superfamily Cerm-13": "Conotoxina_no_clasificada",
    "Superfamily Rmil-02": "Conotoxina_no_clasificada",
    "G119 VD New Superfamily-1": "Conotoxina_no_clasificada",
    "G120 VD New Superfamily-1": "Conotoxina_no_clasificada",
    "G118 VD New Superfamily-1": "Conotoxina_no_clasificada",
    "G127 VD New Superfamily-4": "Conotoxina_no_clasificada",
    "SF1": "Conotoxina_no_clasificada",
    "SF10": "Conotoxina_no_clasificada",
    "I4": "Conotoxina_no_clasificada",
    "I5": "Conotoxina_no_clasificada",

    # exact Ctr values already curated
    "Ctr_2_N": "Conotoxina_no_clasificada",
    "Ctr_6_T (Ctr_P_1)": "Conotoxina_no_clasificada",
    "Ctr_9_T": "Conotoxina_no_clasificada",
    "Ctr_10_TN": "Conotoxina_no_clasificada",
    "Ctr_11_TN (Ctr_S_2)": "Conotoxina_no_clasificada",
    "Ctr_13_TN (Ctr_NL_1)": "Conotoxina_no_clasificada",
    "Ctr_23_T (Ctr_O1_1)": "Conotoxina_no_clasificada",
    "Ctr_24_T": "Conotoxina_no_clasificada",
    "Ctr_26_TN": "Conotoxina_no_clasificada",
    "Ctr_28_N": "Conotoxina_no_clasificada",
    "Ctr_32_T": "Conotoxina_no_clasificada",
    "Ctr_38_TN": "Conotoxina_no_clasificada",
    "Ctr_41_N": "Conotoxina_no_clasificada",
    "Ctr_42_T (Ctr_Di3_1)": "Conotoxina_no_clasificada",
    "Ctr_44_TN": "Conotoxina_no_clasificada",
    "Ctr_48_N (Ctr_M_12)": "Conotoxina_no_clasificada",
    "Ctr_49_N": "Conotoxina_no_clasificada",
    "Ctr_52_TN": "Conotoxina_no_clasificada",
    "Ctr_53_N (Ctr_L_2)": "Conotoxina_no_clasificada",
    "Ctr_54_T": "Conotoxina_no_clasificada",
    "Ctr_55_TN (Ctr_M_1)": "Conotoxina_no_clasificada",
    "Ctr_65_TN": "Conotoxina_no_clasificada",
    "Ctr_68_N": "Conotoxina_no_clasificada",
    "Ctr_69_N (Ctr_I2_2)": "Conotoxina_no_clasificada",
    "Ctr_71_N": "Conotoxina_no_clasificada",
    "Ctr_72_TN (Ctr_Di4_1)": "Conotoxina_no_clasificada",
    "Ctr_75_N": "Conotoxina_no_clasificada",
    "Ctr_77_N": "Conotoxina_no_clasificada",
    "Ctr_82_TN (Ctr_I1_1)": "Conotoxina_no_clasificada",
    "Ctr_84_T": "Conotoxina_no_clasificada",
    "Ctr_85_TN (Ctr_Ikot_3)": "Conotoxina_no_clasificada",
    "Ctr_87_T (Ctr_M_11)": "Conotoxina_no_clasificada",
    "Ctr_92_T (Ctr_O2_1)": "Conotoxina_no_clasificada",
    "Ctr_95_T": "Conotoxina_no_clasificada",
    "Ctr_100_T": "Conotoxina_no_clasificada",
    "Ctr_101_T": "Conotoxina_no_clasificada",
    "Ctr_104_T (Ctr_B1_1)": "Conotoxina_no_clasificada",
    "Ctr_112_T": "Conotoxina_no_clasificada",
    "Ctr_116_T (Ctr_Dip_2)": "Conotoxina_no_clasificada",
    "Ctr_118_T": "Conotoxina_no_clasificada",
    "Ctr_119_T": "Conotoxina_no_clasificada",
    "Ctr_122_T (Ctr_Pre_1)": "Conotoxina_no_clasificada",
    "Ctr_126_N (Ctr_H_2)": "Conotoxina_no_clasificada",
    "Ctr_131_T": "Conotoxina_no_clasificada",
    "Ctr_132_T": "Conotoxina_no_clasificada",
    "Ctr_134_N": "Conotoxina_no_clasificada",
    "Ctr_137_N (Ctr_Ikot_7)": "Conotoxina_no_clasificada",
    "Ctr_143_TN (Ctr_Y2_1) (Insulin)": "Conotoxina_no_clasificada",
    "Ctr_145_TN (Ctr_W_1)": "Conotoxina_no_clasificada",
    "Ctr_147_T": "Conotoxina_no_clasificada",
    "Ctr_148_N": "Conotoxina_no_clasificada",
    "Ctr_149_TN (Ctr_I1_2)": "Conotoxina_no_clasificada",
    "Ctr_153_T": "Conotoxina_no_clasificada",
    "Ctr_156_T (Ctr_Di2_1)": "Conotoxina_no_clasificada",
    "Ctr_157_T": "Conotoxina_no_clasificada",
    "Ctr_160_N": "Conotoxina_no_clasificada",
    "Ctr_D_1": "Conotoxina_no_clasificada",
    "Ctr_Di1_3": "Conotoxina_no_clasificada",
    "Ctr_Di2_2": "Conotoxina_no_clasificada",
    "Ctr_Dip_1": "Conotoxina_no_clasificada",
    "Ctr_Dip_3": "Conotoxina_no_clasificada",
    "Ctr_I3_1": "Conotoxina_no_clasificada",
    "Ctr_Ikot_2": "Conotoxina_no_clasificada",
    "Ctr_Ikot_12": "Conotoxina_no_clasificada",
    "Ctr_K_1": "Conotoxina_no_clasificada",
    "Ctr_M_2": "Conotoxina_no_clasificada",
    "Ctr_M_7": "Conotoxina_no_clasificada",
    "Ctr_N_3": "Conotoxina_no_clasificada",
    "Ctr_O1_3": "Conotoxina_no_clasificada",
    "Ctr_O3_2": "Conotoxina_no_clasificada",
    "Ctr_SF1_1": "Conotoxina_no_clasificada",
    "Ctr_SF2_1": "Conotoxina_no_clasificada",
    "Ctr_SF3_1": "Conotoxina_no_clasificada",
    "Ctr_SF4_1": "Conotoxina_no_clasificada",
    "Ctr_U_1": "Conotoxina_no_clasificada",
    "Ctr_Y_1": "Conotoxina_no_clasificada",

    # Di / Eb / Im / im / Mi / XV exacts
    "Di6": "Conotoxina_no_clasificada",
    "Di6.1": "Conotoxina_no_clasificada",
    "Di6.5": "Conotoxina_no_clasificada",
    "Di6.6": "Conotoxina_no_clasificada",
    "Di6.8": "Conotoxina_no_clasificada",
    "Di6.11": "Conotoxina_no_clasificada",
    "Di6.13": "Conotoxina_no_clasificada",
    "Di7.3": "Conotoxina_no_clasificada",
    "Di7.5": "Conotoxina_no_clasificada",
    "Di14.1": "Conotoxina_no_clasificada",
    "Eb3": "Conotoxina_no_clasificada",
    "Eb3.1": "Conotoxina_no_clasificada",
    "Eb3.2": "Conotoxina_no_clasificada",
    "Eb3.5": "Conotoxina_no_clasificada",
    "Eb5.4": "Conotoxina_no_clasificada",
    "Eb6.6": "Conotoxina_no_clasificada",
    "Eb6.8": "Conotoxina_no_clasificada",
    "Eb14.7": "Conotoxina_no_clasificada",
    "Im11.13": "Conotoxina_no_clasificada",
    "Im11.14": "Conotoxina_no_clasificada",
    "Im14": "Conotoxina_no_clasificada",
    "Im14.1": "Conotoxina_no_clasificada",
    "Im6.1": "Conotoxina_no_clasificada",
    "Im6.4": "Conotoxina_no_clasificada",
    "Im9.12": "Conotoxina_no_clasificada",
    "im001": "Conotoxina_no_clasificada",
    "im002": "Conotoxina_no_clasificada",
    "im017": "Conotoxina_no_clasificada",
    "im020": "Conotoxina_no_clasificada",
    "im021": "Conotoxina_no_clasificada",
    "im031": "Conotoxina_no_clasificada",

    "ViKr111P": "Conotoxina_no_clasificada",
    "ViKr35P": "Conotoxina_no_clasificada",
    "Vi6.7": "Conotoxina_no_clasificada",
    "Vi6.9": "Conotoxina_no_clasificada",
    "Mr3.9": "Conotoxina_no_clasificada",
    "Mr3.10": "Conotoxina_no_clasificada",
    "Mr1.7": "Conotoxina_no_clasificada",
    "lt7b": "Conotoxina_no_clasificada",
    "Mu8.1 (Conus mucronatus)": "Conotoxina_no_clasificada",
    "Vc6.28 prepropeptide": "Conotoxina_no_clasificada",
    "qc23a": "Conotoxina_no_clasificada",
    "MgJ42P": "Conotoxina_no_clasificada",
    "MgJr94P": "Conotoxina_no_clasificada",
    "VeG51P": "Conotoxina_no_clasificada",
    "VeG52P": "Conotoxina_no_clasificada",
    "Bu3.9 protein": "Conotoxina_no_clasificada",
    "Lp7.2": "Conotoxina_no_clasificada",
    "Lt3.6": "Conotoxina_no_clasificada",
    "PIIIH": "Conotoxina_no_clasificada",
    "PIIIG": "Conotoxina_no_clasificada",

    "XV Ca15a": "Conotoxina_no_clasificada",
    "XV Rt15a": "Conotoxina_no_clasificada",
    "XV Rt15b": "Conotoxina_no_clasificada",
    "XV Rt15c": "Conotoxina_no_clasificada",
    "XV S15a": "Conotoxina_no_clasificada",
    "XV Tx15a": "Conotoxina_no_clasificada",
    "XV Vr15a": "Conotoxina_no_clasificada",
    "XV Vr15b": "Conotoxina_no_clasificada",
    "XV Vt15a": "Conotoxina_no_clasificada",
    "XV Vt15b": "Conotoxina_no_clasificada",

    # Mi series
    **{f"Mi{i:03d}": "Conotoxina_no_clasificada" for i in range(1, 49) if i not in {40, 42, 45, 46}},
}

# v12 overlay from user-provided homologation table
EXACT_MAP.update({
    'superfamily A': 'A',
    'Superfamily A': 'A',
    'superfamily A2': 'A2',
    'superfamily B1': 'B1',
    'superfamily B2': 'B2',
    'superfamily C': 'C',
    'superfamily con-ikot-ikot': 'con-ikot-ikot',
    'superfamily conkunitzin 3': 'Conkunitzin',
    'superfamily conkunitzin-7': 'Conkunitzin',
    'superfamily conkunitzin-3': 'Conkunitzin',
    'superfamily conkunitzin-4': 'Conkunitzin',
    'superfamily conkunitzin 1': 'Conkunitzin',
    'superfamily conkunitzin-6': 'Conkunitzin',
    'Superfamily Ggeo03-2': 'Conotoxina_no_clasificada',
    'Superfamily Cerm-07': 'Conotoxina_no_clasificada',
    'Superfamily SF1': 'Conotoxina_no_clasificada',
    'Superfamily Cerm-11': 'Conotoxina_no_clasificada',
    'Superfamily SF10': 'Conotoxina_no_clasificada',
    'Superfamily Cerm-01': 'Conotoxina_no_clasificada',
    'Superfamily Cerm-06': 'Conotoxina_no_clasificada',
    'Superfamily Cerm-10': 'Conotoxina_no_clasificada',
    'Superfamily Rimp-05': 'Conotoxina_no_clasificada',
    'Superfamily Ggeo03-1': 'Conotoxina_no_clasificada',
    'Cysteine-rich protein 1': 'Conotoxina_no_clasificada',
    'Secreted protein': 'Conotoxina_no_clasificada',
    'superfamily D': 'D',
    'superfamily E': 'E',
    'superfamily F': 'F',
    'superfamily H': 'H',
    'superfamily I1': 'I1',
    'Superfamily I1': 'I1',
    'superfamily I2': 'I2',
    'Superfamily I2': 'I2',
    'superfamily I3': 'I3',
    'superfamily I4': 'I4',
    'superfamily I5': 'I5',
    'superfamily J': 'J',
    'superfamily K': 'K',
    'superfamily L': 'L',
    'superfamily M': 'M',
    'Superfamily M': 'M',
    'Neuropeptide F': 'NO_conotoxin',
    'Hyaluronidase (EC 3.2.1.35)': 'NO_conotoxin',
    'peptidylglycine monooxygenase (EC 1.14.17.3)': 'NO_conotoxin',
    'Thyrostimulin alpha': 'NO_conotoxin',
    'Thyrostimulin beta': 'NO_conotoxin',
    'Lysozyme': 'NO_conotoxin',
    'Thioredoxin domain-containing protein 5': 'NO_conotoxin',
    'superfamily O1': 'O1',
    'Superfamily O1': 'O1',
    'superfamily O2': 'O2',
    'Superfamily O2': 'O2',
    'superfamily O3': 'O3',
    'superfamily P': 'P',
    'superfamily Q': 'Q',
    'superfamily S': 'S',
    'Superfamily S': 'S',
    'superfamily T': 'T',
    'superfamily V': 'V',
    'superfamily W': 'W',
    'superfamily Z': 'Z',
    'Superfamily Z': 'Z',
})
VALID_FINAL.update({'A2', 'I4', 'I5'})

EXACT_MAP_CANON = {canon_key(k): v for k, v in EXACT_MAP.items()}

RE_CONKUNITZIN_VARIANT = re.compile(r'^conkunitzin[\s-]*\d+$', re.I)
RE_CERM_VARIANT = re.compile(r'^Cerm-\d+$', re.I)
RE_RIMP_VARIANT = re.compile(r'^Rimp-\d+$', re.I)
RE_GGEO_VARIANT = re.compile(r'^Ggeo\d+-\d+$', re.I)
RE_SF_NUMBER = re.compile(r'^SF\d+$', re.I)
RE_I_NUMBER_UNCLASSIFIED = re.compile(r'^I([4-9]\d*|\d{2,})$', re.I)
RE_MI_SERIES = re.compile(r'^Mi\d+$', re.I)
RE_IM_SERIES = re.compile(r'^(?:im\d+|Im\d+(?:\.\d+)?)$', re.I)
RE_DI_SERIES = re.compile(r'^Di\d+(?:\.\d+)?$', re.I)
RE_EB_SERIES = re.compile(r'^Eb\d+(?:\.\d+)?$', re.I)
RE_XV_SERIES = re.compile(r'^XV\s+\S.*$', re.I)

RE_VD_RS_SUPERFAMILY = re.compile(rf'^G\d+\s+(?:VD|RS)\s+Superfamily\s+({FAMILY_ALT})\b', re.I)
RE_LEADING_CODE_SUPERFAMILY = re.compile(rf'^({FAMILY_ALT})[\s-]+superfamily\b', re.I)
RE_LEADING_CODE_UNDERSCORE = re.compile(rf'^({FAMILY_ALT})(?:_|[\s-]*_)', re.I)
RE_CODE_PEPTIDE = re.compile(rf'^({FAMILY_ALT})\s*-\s*peptide\b', re.I)
RE_BARE_CODE_DASH = re.compile(r'^([A-Z][A-Z0-9]*)-\s*$', re.I)  # allows B-
RE_ALPHA = re.compile(r'^Alpha(?:-like)?(?:\b.*)?$', re.I)

NO_CONOTOXIN_PATTERNS = [
    re.compile(r'^Insulin(?:-like|-related)?(?:\s+protein(?:\s+\d+)?)?$', re.I),
    re.compile(r'^Hormone superfamily prohormone-\d+$', re.I),
    re.compile(r'^Protein disulfide[- ]isomerase(?:\b.*)?$', re.I),
    re.compile(r'^protein disulfide[- ]isomerase(?:\b.*)?$', re.I),
    re.compile(r'^Peptidyl-prolyl cis-trans isomerase(?:\b.*)?$', re.I),
    re.compile(r'^Peptidylglycine alpha(?: |-)?amidating monooxygenase(?:\b.*)?$', re.I),
    re.compile(r'^glutaminyl-peptide cyclotransferase(?:\b.*)?$', re.I),
    re.compile(r'^Glutaminyl-peptide cyclotransferase(?:\b.*)?$', re.I),
    re.compile(r'^procollagen-proline 4-dioxygenase(?:\b.*)?$', re.I),
    re.compile(r'^Endoplasmic reticulum oxidoreductin isoform X\d+$', re.I),
    re.compile(r'^Endoplasmic reticulum chaperone BIP(?:\b.*)?$', re.I),
    re.compile(r'^Serine protease inhibitor$', re.I),
    re.compile(r'^Serine protease$', re.I),
    re.compile(r'^Serine/threonine protein kinase \d+$', re.I),
    re.compile(r'^non-specific serine/threonine protein kinase(?:\b.*)?$', re.I),
    re.compile(r'^Potassium voltage-gated channel Kv1\.1-like protein$', re.I),
    re.compile(r'^Vitamin-K dependent carboxylase$', re.I),
    re.compile(r'^vitamin-K-epoxide reductase(?:\b.*)?$', re.I),
    re.compile(r'^Vacuolar proton pump subunit B(?:\b.*)?$', re.I),
    re.compile(r'^Tubulin beta chain$', re.I),
    re.compile(r'^High frequency protein 2$', re.I),
    re.compile(r'^F-ATPase protein 6$', re.I),
    re.compile(r'^Kazal protease inhibitor$', re.I),
    re.compile(r'^Nicotinic acetylcholine receptor alpha 2-like protein$', re.I),
]

def norm_uniprot(x: object) -> str:
    s = normalize_space(x)
    if not s or is_empty_token(s):
        return ""

    ck = canon_key(s)
    if ck in EXACT_MAP_CANON:
        return EXACT_MAP_CANON[ck]

    if ck in {canon_key(v) for v in TO_EMPTY}:
        return ""

    # already final
    if s in VALID_FINAL:
        return s

    # explicit proteins / non-conotoxins
    for pat in NO_CONOTOXIN_PATTERNS:
        if pat.match(s):
            if re.match(r'^Insulin(?:-like|-related)?(?:\s+protein(?:\s+\d+)?)?$', s, re.I):
                return "Insulin"
            return "NO_conotoxin"

    # broad curated families
    if RE_ALPHA.match(s):
        return "A"

    if s.casefold() == "superfamily conantokin":
        return "Conantokin"
    if s.startswith("Conantokin-F peptide"):
        return "Conantokin"

    # names by prefix
    if s.startswith("Conodipine"):
        return "Conodipine"
    if s.startswith("Conomarphin"):
        return "Conomarphin"
    if s.startswith("Contryphan"):
        return "Contryphan"
    if s.startswith("Conorfamide"):
        return "Conorfamide"
    if s == "Conoporin":
        return "Conoporin"
    if s.lower().startswith("con-ikot-ikot"):
        return "con-ikot-ikot"
    if "Conantokin__" in s:
        return "Conantokin"
    if "Conkunitzin__" in s:
        return "Conkunitzin"
    if "Contulakin" in s or s.casefold() == "contulakin":
        return "Contulakin"

    # generic unclassified families
    if RE_CERM_VARIANT.match(s) or RE_RIMP_VARIANT.match(s) or RE_GGEO_VARIANT.match(s) or RE_SF_NUMBER.match(s) or RE_I_NUMBER_UNCLASSIFIED.match(s):
        return "Conotoxina_no_clasificada"

    # broad series
    if s.startswith("Ctr_"):
        return "Conotoxina_no_clasificada"
    if RE_MI_SERIES.match(s) or RE_IM_SERIES.match(s) or RE_DI_SERIES.match(s) or RE_EB_SERIES.match(s) or RE_XV_SERIES.match(s):
        return "Conotoxina_no_clasificada"

    # VD/RS Superfamily X
    m = RE_VD_RS_SUPERFAMILY.match(s)
    if m:
        return m.group(1).upper()

    # X-superfamily / X superfamily
    m = RE_LEADING_CODE_SUPERFAMILY.match(s)
    if m:
        code = m.group(1).upper()
        return "B1" if code == "B-" else code

    # X_....
    m = RE_LEADING_CODE_UNDERSCORE.match(s)
    if m:
        return m.group(1).upper()

    # X- peptide
    m = RE_CODE_PEPTIDE.match(s)
    if m:
        return m.group(1).upper()

    # B- / O2- / p-
    m = RE_BARE_CODE_DASH.match(s)
    if m:
        code = m.group(1).upper()
        if code == "B":
            return "B1"
        if code == "P":
            return "P"
        if code in VALID_FINAL:
            return code

    return s

def normalize_ncbi_value(s: object) -> tuple[str, str]:
    """
    Normaliza NCBI usando únicamente la etiqueta NCBI.

    El catálogo global contiene 29 etiquetas: 28 ya canónicas y una variante
    de capitalización ("conkunitzin").
    """
    t = normalize_space(s)
    if not t or is_empty_token(t):
        return "", "empty"

    if t in VALID_FINAL:
        return t, "already_canonical"

    if canon_key(t) == "conkunitzin":
        return "Conkunitzin", "case_normalization_conkunitzin"

    return t, "unresolved"


# ---------------------------------------------------------------------------
# Normalización global v3
#
# Principio:
# - Solo se asigna una superfamilia cuando el nombre de UniProt la expresa.
# - Los nombres de conotoxinas/conopéptidos sin superfamilia explícita se
#   conservan como evidencia conotoxínica no clasificada.
# - Las proteínas celulares/enzimáticas explícitas se marcan NO_conotoxin.
# - "Uncharacterized protein" se deja vacío porque no aporta evidencia a favor
#   ni en contra de una conotoxina.
# - No se usa ConoServer, NCBI ni Signal para convertir una etiqueta UniProt
#   ambigua en una superfamilia: hacerlo inflaría artificialmente la
#   concordancia entre fuentes.
# ---------------------------------------------------------------------------

NORMALIZATION_VERSION = "global_2026-07-20_v6_project_nomenclature"

FAMILY_CODES_V3 = sorted(
    FINAL_SUPERFAMILIES | {"A2", "I4", "I5"},
    key=lambda value: (-len(value), value),
)
FAMILY_ALT_V3 = "|".join(re.escape(value) for value in FAMILY_CODES_V3)

RE_CONOTOXIN_SUPERFAMILY_V3 = re.compile(
    rf"^Conotoxin\s+superfamily\s+({FAMILY_ALT_V3})\b", re.I
)
RE_CONOTOXIN_CODE_SUPERFAMILY_V3 = re.compile(
    rf"^Conotoxin\s+({FAMILY_ALT_V3})\s+superfamily\b", re.I
)
RE_LEADING_SUPERFAMILY_V3 = re.compile(
    rf"^({FAMILY_ALT_V3})\s+superfamily\b", re.I
)
RE_CODE_CONOTOXIN_V3 = re.compile(
    rf"^({FAMILY_ALT_V3})\s*-\s*conotoxin\b", re.I
)
RE_CODE_CONOTOXIN_SPACE_V3 = re.compile(
    rf"^({FAMILY_ALT_V3})\s+conotoxin\b", re.I
)
RE_CODE_PREPROPEPTIDE_V3 = re.compile(
    rf"^({FAMILY_ALT_V3})_[A-Za-z0-9].*\b(?:prepropeptide|precursor)\b",
    re.I,
)
RE_VD_RS_SUPERFAMILY_V3 = re.compile(
    rf"^G\d+[\s_]+(?:VD|RS)[\s_]+(?:New[\s_]+)?"
    rf"Superfamily[\s_]+({FAMILY_ALT_V3})(?:-variant\d+)?\b",
    re.I,
)
RE_EMBEDDED_SUPERFAMILY_V3 = re.compile(
    rf"\b({FAMILY_ALT_V3})-superfamily\s+protein\b", re.I
)

NON_CONOTOXIN_PATTERNS_V3 = [
    re.compile(r"^procollagen-proline 4-dioxygenase\b", re.I),
    re.compile(r"^Endoplasmic reticulum chaperone BIP\b", re.I),
    re.compile(r"^(?:protein )?disulfide-isomerase\b", re.I),
    re.compile(r"^Peptidyl-prolyl cis-trans isomerase\b", re.I),
    re.compile(r"^Metalloendopeptidase\b", re.I),
    re.compile(r"^[Gg]lutaminyl-peptide cyclotransferase\b", re.I),
    re.compile(r"^Peptidylglycine(?: alpha amidating| alpha-amidating| )?"
               r"monooxygenase\b", re.I),
    re.compile(r"^peptidylglycine monooxygenase\b", re.I),
    re.compile(r"^Kazal protease inhibitor\b", re.I),
    re.compile(r"^Serine protease(?: inhibitor)?\b", re.I),
    re.compile(r"^non-specific serine/threonine protein kinase\b", re.I),
    re.compile(r"^Nicotinic acetylcholine receptor\b", re.I),
    re.compile(r"^Potassium voltage-gated channel\b", re.I),
    re.compile(r"^Hyaluronidase\b", re.I),
    re.compile(r"^Thioredoxin domain-containing protein\b", re.I),
    re.compile(r"^(?:78 kDa|78Kda) glucose[- ]regulated protein\b", re.I),
    re.compile(r"^Cytochrome\b", re.I),
    re.compile(r"^NADH-ubiquinone oxidoreductase\b", re.I),
    re.compile(r"^F-ATPase protein\b", re.I),
    re.compile(r"^Vacuolar proton pump\b", re.I),
    re.compile(r"^Vitamin-K dependent carboxylase\b", re.I),
    re.compile(r"^vitamin-K-epoxide reductase\b", re.I),
    re.compile(r"^Elongation factor\b", re.I),
    re.compile(r"^Ferritin\b", re.I),
    re.compile(r"^Arginine kinase\b", re.I),
    re.compile(r"^Lysozyme\b", re.I),
    re.compile(r"^Neuropeptide F$", re.I),
    re.compile(r"^Neuropeptide prohormone\b", re.I),
    re.compile(r"^Hormone superfamily prohormone\b", re.I),
    re.compile(r"^Thyrostimulin\b", re.I),
    re.compile(r"^thiol oxidase\b", re.I),
    re.compile(r"^High frequency protein 2$", re.I),
    re.compile(r"^Cysteine-rich venom protein\b", re.I),
]

RE_GENERIC_CONOTOXIN_V3 = re.compile(
    r"^(?:Conotoxin|Putative conotoxin|Four-loop conotoxin|"
    r"Macroconotoxin|Conopeptide|Linear conopeptide)\b",
    re.I,
)
RE_PHARMACOLOGICAL_CONOTOXIN_V3 = re.compile(
    r"(?:^|[- ])(?:alpha|beta|gamma|delta|epsilon|omega|mu|kappa|"
    r"rho|sigma|tau|chi)-?conotoxin",
    re.I,
)
RE_NAMED_CONOPEPTIDE_V3 = re.compile(
    r"^(?:"
    r"Mr(?:\d+\.\d+|_precursor_\d+|1\.1_precursor)|"
    r"Mi\d+|Im\d+(?:\.\d+)?|Di\d+(?:\.\d+)?|Eb\d+(?:\.\d+)?|"
    r"Tx\d+(?:\.\d+)?|ViKr\d+P|MgJr?\d+P|VeG\d+P|p21b|"
    r"PIII[GH]|Bu\d+\.\d+\s+protein|TxMMSK-\d+"
    r")\b",
    re.I,
)


# IMPORTANTE:
# Las clases peptídicas de UniProt (por ejemplo, Contryphan) se conservan
# literalmente en la columna de origen. Su correspondencia con una
# superfamilia génica se aplica únicamente durante la clasificación.
#
def normalize_uniprot_v3(x: object) -> tuple[str, str]:
    """Return (normalized_label, applied_rule)."""
    s = normalize_space(x)
    if not s or is_empty_token(s):
        return "", "empty"

    # Corregir antes de las reglas históricas: "E1" en "Conantokin E1"
    # es parte del nombre del péptido y no una superfamilia génica.
    # También se normalizan las variantes textuales conantokin-like.
    if re.search(r"\bconantokin(?:[-\s]?like)?\b", s, re.I):
        return "Conantokin", "special_conantokin_early"

    # Preserve all prior curated decisions that already produce a valid value.
    previous = norm_uniprot(s)
    if previous in VALID_FINAL or previous == "":
        return previous, "existing_curated_rule"

    # New named venom-peptide classes present in the updated UniProt dataset.
    if re.search(r"\bconsomatin\b", s, re.I):
        return "Consomatin", "special_consomatin"
    if re.match(r"^ConoCAP\b", s, re.I):
        return "ConoCAP", "special_conocap"
    if re.match(r"^Elevenin\b", s, re.I):
        return "Elevenin", "special_elevenin"
    if re.match(r"^Scratcher peptide\b", s, re.I):
        return "Scratcher", "special_scratcher"

    # Existing special classes, expanded to cover wording variants.
    if re.search(r"contryphan", s, re.I):
        return "Contryphan", "special_contryphan"
    if re.search(r"conopressin|conophysin", s, re.I):
        return "Conopressin", "special_conopressin"
    if re.search(r"conkunitzin", s, re.I):
        return "Conkunitzin", "special_conkunitzin"
    if re.search(r"conantokin", s, re.I):
        return "Conantokin", "special_conantokin"
    if re.search(r"con-ikot-ikot", s, re.I):
        return "con-ikot-ikot", "special_con_ikot_ikot"
    if re.search(r"conodipine", s, re.I):
        return "Conodipine", "special_conodipine"
    if re.search(r"conomarphin", s, re.I):
        return "Conomarphin", "special_conomarphin"
    if re.search(r"conorfamide", s, re.I):
        return "Conorfamide", "special_conorfamide"
    if re.match(r"^Conoporin\b", s, re.I):
        return "Conoporin", "special_conoporin"
    if re.match(r"^Con-Ins\b", s, re.I):
        return "Insulin", "special_insulin"
    if re.match(r"^Insulin(?:-like|-related)?\b", s, re.I):
        return "Insulin", "special_insulin"

    # Superfamily explicitly encoded in the protein name.
    embedded = RE_EMBEDDED_SUPERFAMILY_V3.search(s)
    if embedded:
        return embedded.group(1).upper(), "embedded_superfamily"

    explicit_family_patterns = [
        (RE_CONOTOXIN_SUPERFAMILY_V3, "conotoxin_superfamily"),
        (RE_CONOTOXIN_CODE_SUPERFAMILY_V3, "conotoxin_code_superfamily"),
        (RE_LEADING_SUPERFAMILY_V3, "leading_superfamily"),
        (RE_CODE_CONOTOXIN_V3, "code_conotoxin"),
        (RE_CODE_CONOTOXIN_SPACE_V3, "code_conotoxin"),
        (RE_CODE_PREPROPEPTIDE_V3, "code_prepropeptide"),
        (RE_VD_RS_SUPERFAMILY_V3, "vd_rs_superfamily"),
    ]
    for pattern, rule in explicit_family_patterns:
        match = pattern.match(s)
        if match:
            return match.group(1).upper(), rule

    # Explicitly conotoxin-like, but without an assignable established family.
    if re.search(r"Conotoxin-like\s+unassigned\s+superfamily", s, re.I):
        return "Conotoxina_no_clasificada", "unassigned_superfamily"
    if re.search(r"Conotoxin\s+unclassified\s+superfamily", s, re.I):
        return "Conotoxina_no_clasificada", "unclassified_superfamily"
    if re.search(r"New\s+Superfamily-\d+", s, re.I):
        return "Conotoxina_no_clasificada", "new_superfamily"
    if re.match(
        r"^(?:Superfamily|Conotoxin Superfamily)\s+"
        r"(?:SF\d+|Cerm-\d+|Rimp-\d+|Tand-\d+|Ggeo\d+-\d+|"
        r"Pmag-\d+|Rmil-\d+)\b",
        s,
        re.I,
    ):
        return "Conotoxina_no_clasificada", "named_unclassified_superfamily"
    if re.match(r"^Superfamily conotoxin$", s, re.I):
        return "Conotoxina_no_clasificada", "generic_conotoxin"

    # Uninformative label: do not treat it as evidence for either class.
    if canon_key(s) == canon_key("Uncharacterized protein"):
        return "", "uninformative_uncharacterized"

    # Explicit proteins or enzymes not annotated as conotoxins.
    for pattern in NON_CONOTOXIN_PATTERNS_V3:
        if pattern.match(s):
            return "NO_conotoxin", "explicit_non_conotoxin"

    # Pharmacological names do not by themselves define a gene superfamily.
    if RE_PHARMACOLOGICAL_CONOTOXIN_V3.search(s):
        return "Conotoxina_no_clasificada", "pharmacological_conotoxin"

    # Named conotoxin/conopeptide without an explicit gene-superfamily label.
    if RE_GENERIC_CONOTOXIN_V3.match(s):
        return "Conotoxina_no_clasificada", "generic_or_named_conotoxin"
    if re.match(
        r"^G\d+[\s_]+(?:VD|RS)[\s_]+.*(?:conopeptide|precursor)",
        s,
        re.I,
    ):
        return "Conotoxina_no_clasificada", "unnamed_conopeptide"
    if RE_NAMED_CONOPEPTIDE_V3.match(s):
        return "Conotoxina_no_clasificada", "named_conopeptide"
    if s in {"Cysteine-rich protein 1", "Secreted protein"}:
        return "Conotoxina_no_clasificada", "generic_conotoxin_like"

    return s, "unresolved"


def norm_uniprot_v3(x: object) -> str:
    return normalize_uniprot_v3(x)[0]
