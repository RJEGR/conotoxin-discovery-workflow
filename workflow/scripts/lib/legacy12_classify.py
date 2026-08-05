#!/usr/bin/env python3
# =============================================================================
# legacy12_classify.py  -  Stage 7 confidence classification (vendored)
# =============================================================================
# VENDORED VERBATIM from the legacy Stage 7 script:
#     context/Legacy_stage7/12_clasificar_por_certidumbre_v7.py
#
# Everything below the header is an unmodified copy of that script's
# confidence-classification layer (from the top of the file down to, but excluding,
# `def reorder_columns`, which belongs to the legacy I/O layer that
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

import argparse
import csv
from pathlib import Path
import re
from typing import Dict, List, Tuple
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


EMPTY_TOKENS = {"", "nan", "none", "null", "na", "n/a", "-", "—"}


def _norm_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    value = str(value).strip()
    return re.sub(r"\s+", " ", value)


def _normalize_superfamily_label(value: object) -> str:
    """Normaliza nomenclatura para comparación interna."""
    value = _norm_text(value)
    if not value or value.lower() in EMPTY_TOKENS:
        return ""

    if value.upper() == "NO_CONOTOXIN":
        return "NO_CONOTOXIN"

    # Para comparación, tratar esto como ausencia de clasificación
    if value.upper() == "CONOTOXINA_NO_CLASIFICADA":
        return ""

    # Casos específicos antes de convertir a mayúsculas
    if value == "MSTLGMTLL-":
        return value

    value_up = value.upper()

    patterns = [
        (r'^(.+?)\s*-?\s*SUPERFAMILY\s*$', r'\1'),
        (r'^SUPERFAMILY\s+(.+?)$', r'\1'),
        (r'^(.+?)\s*-?\s*CONOTOXIN\s*$', r'\1'),
        (r'^(.+?)\s*-?\s*PEPTIDE\s*$', r'\1'),
        (r'^ALPHA(?:-LIKE)?(?:\b.*)?$', 'A'),
    ]

    for pattern, replacement in patterns:
        if re.match(pattern, value_up, re.IGNORECASE):
            value_up = re.sub(pattern, replacement, value_up, flags=re.IGNORECASE).strip()
            break

    equivalencias = {
        "INSULIN": "INSULIN",
        "INSULIN-RELATED PROTEIN": "INSULIN",
        "INSULIN-LIKE PROTEIN": "INSULIN",
        "INSULIN-LIKE PROTEIN 1": "INSULIN",
        "CONKUNITZIN": "CONKUNITZIN",
        "CONANTOKIN": "CONANTOKIN",
    }

    return equivalencias.get(value_up, value_up)


PEPTIDE_CLASS_TO_GENE_SUPERFAMILY = {
    "CONTRYPHAN": "O2",
    "CONANTOKIN": "B1",
    "CONOMARPHIN": "M",
    "CONSOMATIN": "C",
    "SCRATCHER": "J",
}

PROJECT_DIRECT_EQUIVALENCES = {
    "Y2": "INSULIN",
    "Y3": "Y",
    "H2": "H",
    "N2": "N",
    "A2": "A",
}


def normalize_superfamily_for_compare_with_note(
    value: object,
) -> Tuple[str, str]:
    """
    Devuelve (etiqueta_comparable, nota_de_mapeo).

    Las columnas originales no se modifican. Los mapeos se aplican
    únicamente para calcular compatibilidad y consenso.
    """
    normalized = _normalize_superfamily_label(value)

    if not normalized:
        return "", ""

    # Conantokin E1 y conantokin-like se interpretan como la clase
    # Conantokin; no se aplica una regla global E1 -> B1.
    if (
        normalized == "CONANTOKIN"
        or normalized.startswith("CONANTOKIN E1")
        or normalized.startswith("CONANTOKIN-E1")
        or normalized.startswith("CONANTOKIN-LIKE")
        or normalized.startswith("CONANTOKIN LIKE")
    ):
        return "B1", f"{normalized}_a_B1"

    target = PEPTIDE_CLASS_TO_GENE_SUPERFAMILY.get(normalized)
    if target:
        return target, f"{normalized}_a_{target}"

    target = PROJECT_DIRECT_EQUIVALENCES.get(normalized)
    if target:
        return target, f"{normalized}_a_{target}"

    return normalized, ""


def normalize_superfamily_for_compare(value: object) -> str:
    """Devuelve solo la etiqueta comparable."""
    normalized, _ = normalize_superfamily_for_compare_with_note(value)
    return normalized


def normalize_superfamily_for_output(value: object) -> str:
    """Normaliza nomenclatura para mostrar en el Excel de salida."""
    value = _norm_text(value)
    if not value or value.lower() in EMPTY_TOKENS:
        return ""

    if value.upper() == "NO_CONOTOXIN":
        return "NO_conotoxin"

    if value.upper() == "CONOTOXINA_NO_CLASIFICADA":
        return "Conotoxina_no_clasificada"

    if value == "MSTLGMTLL-":
        return value

    # Mantener algunas capitalizaciones informativas
    exact_keep = {
        "Conopressin", "Conantokin", "Conkunitzin", "Contulakin",
        "Conodipine", "Conomarphin", "Contryphan", "Conoporin",
        "Conorfamide", "con-ikot-ikot", "NO_conotoxin", "Conotoxina_no_clasificada"
    }
    if value in exact_keep:
        return value

    # Limpiezas suaves
    patterns = [
        (r'^(.+?)\s*-?\s*superfamily\s*$', r'\1'),
        (r'^superfamily\s+(.+?)$', r'\1'),
        (r'^(.+?)\s*-?\s*conotoxin\s*$', r'\1'),
        (r'^(.+?)\s*-?\s*peptide\s*$', r'\1'),
    ]
    out = value
    for pattern, replacement in patterns:
        if re.match(pattern, out, re.IGNORECASE):
            out = re.sub(pattern, replacement, out, flags=re.IGNORECASE).strip()
            break

    if re.match(r'^alpha(?:-like)?(?:\b.*)?$', out, re.IGNORECASE):
        return "A"

    # Normalizaciones de términos frecuentes
    lower = out.lower()
    if lower in {"insulin", "insulin-related protein", "insulin-like protein", "insulin-like protein 1"}:
        return "Insulin"
    if lower == "conkunitzin":
        return "Conkunitzin"
    if lower == "conantokin":
        return "Conantokin"

    # Si es un código de superfamilia, devolverlo en mayúsculas
    if re.match(r'^(?:A|B|B1|B2|B3|C|D|E|E1|F|G|G2|H|I|I1|I2|I3|J|K|L|M|N|O|O1|O2|O3|P|Q|R|S|T|U|V|W|Y|Z)$', out, re.IGNORECASE):
        return out.upper()

    return out


def resolve_compatible_superfamilies(values: List[str]) -> Tuple[str, str]:
    """
    Resuelve el consenso nomenclatural sin modificar las fuentes.

    Reglas:
    - B, O e I son comodines de una única subfamilia específica.
    - I2 es compatible como etiqueta legacy con I4 o I5.
    - M2 + M se resuelve como M, con certidumbre máxima MEDIA.
    - O4 adopta una única familia O/O1/O2/O3, con máximo MEDIA.
    - Y2/Y3/H2/N2/A2 se transforman internamente mediante la tabla
      PROJECT_DIRECT_EQUIVALENCES.
    """
    clean = [normalize_superfamily_for_compare(value) for value in values]
    clean = [
        value for value in clean
        if value not in {"", "NO_CONOTOXIN"}
    ]

    if not clean:
        return "", "sin_anotaciones"

    unique = set(clean)

    if len(unique) == 1:
        return clean[0], "coincidencia_exacta"

    # M2 es compatible únicamente con M.
    if "M2" in unique:
        if unique.issubset({"M2", "M"}):
            return "M", "compatibilidad_condicional_M2_M"
        return "", "conflicto_M2_con_otra_familia"

    # O4 funciona como comodín condicional dentro del grupo O.
    if "O4" in unique:
        allowed = {"O4", "O", "O1", "O2", "O3"}
        if not unique.issubset(allowed):
            return "", "conflicto_O4_con_familia_no_O"

        specifics = unique.intersection({"O1", "O2", "O3"})
        if len(specifics) > 1:
            return "", "conflicto_O4_multiples_O_especificas"
        if len(specifics) == 1:
            return next(iter(specifics)), "compatibilidad_condicional_O4_O"
        if "O" in unique:
            return "O", "compatibilidad_condicional_O4_O"
        return "O4", "coincidencia_exacta"

    # B genérico.
    b_group = {"B", "B1", "B2", "B3"}
    if unique.issubset(b_group):
        specifics = unique - {"B"}
        if len(specifics) == 1:
            return next(iter(specifics)), "compatibilidad_B_generica"
        return "", "conflicto_subfamilias_B"

    # O genérico.
    o_group = {"O", "O1", "O2", "O3"}
    if unique.issubset(o_group):
        specifics = unique - {"O"}
        if len(specifics) == 1:
            return next(iter(specifics)), "compatibilidad_O_generica"
        return "", "conflicto_subfamilias_O"

    # I genérico e I2 legacy.
    i_group = {"I", "I1", "I2", "I3", "I4", "I5"}
    if unique.issubset(i_group):
        if "I4" in unique:
            if unique.issubset({"I", "I2", "I4"}):
                return "I4", "compatibilidad_legacy_I2_I4"
            return "", "conflicto_subfamilias_I"

        if "I5" in unique:
            if unique.issubset({"I", "I2", "I5"}):
                return "I5", "compatibilidad_legacy_I2_I5"
            return "", "conflicto_subfamilias_I"

        specifics = unique - {"I"}
        if len(specifics) == 1:
            return next(iter(specifics)), "compatibilidad_I_generica"

        return "", "conflicto_subfamilias_I"

    return "", "conflicto_entre_superfamilias"


def compatible_pair(value_a: str, value_b: str) -> Tuple[bool, str, str]:
    """Devuelve (compatibles, consenso, relación) para dos etiquetas."""
    consensus, relation = resolve_compatible_superfamilies(
        [value_a, value_b]
    )
    return bool(consensus), consensus, relation


def relation_suffix(relation: str) -> str:
    mapping = {
        "coincidencia_exacta": "coinciden",
        "compatibilidad_B_generica": "compatibilidad_B",
        "compatibilidad_O_generica": "compatibilidad_O",
        "compatibilidad_I_generica": "compatibilidad_I",
        "compatibilidad_legacy_I2_I4": "compatibilidad_legacy_I2_I4",
        "compatibilidad_legacy_I2_I5": "compatibilidad_legacy_I2_I5",
        "compatibilidad_condicional_M2_M": "compatibilidad_condicional_M2_M",
        "compatibilidad_condicional_O4_O": "compatibilidad_condicional_O4_O",
    }
    return mapping.get(relation, relation)


CONDITIONAL_MEDIA_RELATIONS = {
    "compatibilidad_condicional_M2_M",
    "compatibilidad_condicional_O4_O",
}

HARD_PROJECT_CONFLICT_RELATIONS = {
    "conflicto_M2_con_otra_familia",
    "conflicto_O4_con_familia_no_O",
    "conflicto_O4_multiples_O_especificas",
}


def extract_sample_name(file_path: Path) -> str:
    name = file_path.stem
    for suffix in ["_filtrado_signalp_reacomodo_positivo_superfamilia",
                   "_filtrado_signalp_reacomodo_negativo_superfamilia"]:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def _classify_canonical_sequence(
    is_positive: bool,
    conosorter: str,
    conoserver: str,
    ncbi: str,
    uniprot: str,
) -> Tuple[str, str, str]:
    """
    Aplica las reglas de clasificación y retorna:
        (certidumbre, regla, consenso)

    Las cuatro etiquetas se comparan sin reescribir las columnas originales.
    Las compatibilidades B/I solo se aplican al cálculo del consenso.
    """
    conosorter = normalize_superfamily_for_compare(conosorter)
    conoserver = normalize_superfamily_for_compare(conoserver)
    ncbi = normalize_superfamily_for_compare(ncbi)
    uniprot = normalize_superfamily_for_compare(uniprot)

    # Regla histórica específica del proyecto.
    if (
        conosorter == "N"
        and conoserver == "MSTLGMTLL-"
        and ncbi in {"", "N"}
        and uniprot in {"", "N", "MSTLGMTLL-"}
    ):
        return (
            "MEDIA",
            "Especial_N_MSTLGMTLL_apoyo_parcial",
            "MSTLGMTLL-",
        )

    db_values = [conoserver, ncbi, uniprot]
    db_non_empty = [value for value in db_values if value != ""]

    no_conotoxin_count = sum(
        value == "NO_CONOTOXIN" for value in db_non_empty
    )
    superfamily_annotations = [
        value for value in db_non_empty
        if value != "NO_CONOTOXIN"
    ]

    if no_conotoxin_count > 0:
        if not superfamily_annotations:
            return (
                "SIN_CLASIFICAR",
                "Solo_NO_conotoxin",
                "NO_CONOTOXIN",
            )
        return (
            "BAJA",
            "NO_conotoxin_conflicto_con_superfamilias",
            "",
        )

    db_consensus, db_relation = resolve_compatible_superfamilies(
        superfamily_annotations
    )

    if is_positive:
        if not superfamily_annotations:
            return "BAJA", "Solo_ConoSorter", conosorter

        pair_results = [
            compatible_pair(conosorter, value)
            for value in superfamily_annotations
        ]
        compatible_pairs = [
            result for result in pair_results if result[0]
        ]
        n_compatible = len(compatible_pairs)

        all_consensus, all_relation = resolve_compatible_superfamilies(
            [conosorter, *superfamily_annotations]
        )

        # Las reglas M2/O4 no pueden rescatarse mediante una coincidencia
        # parcial cuando existe una familia externa incompatible.
        if all_relation in HARD_PROJECT_CONFLICT_RELATIONS:
            return "BAJA", all_relation, ""

        # M2/M y O4/O son compatibilidades condicionadas: aunque todas las
        # fuentes sean compatibles, la certidumbre máxima es MEDIA.
        if (
            all_consensus
            and all_relation in CONDITIONAL_MEDIA_RELATIONS
        ):
            return (
                "MEDIA",
                relation_suffix(all_relation) + "_max_MEDIA",
                all_consensus,
            )

        # Todas las anotaciones presentes son compatibles entre sí.
        if all_consensus:
            suffix = relation_suffix(all_relation)

            if len(superfamily_annotations) == 1:
                if all_relation == "coincidencia_exacta":
                    rule = "ConoSorter+1DB_coinciden_otras_vacias"
                else:
                    rule = (
                        "ConoSorter+1DB_compatible_otras_vacias_"
                        + suffix
                    )
                return "ALTA", rule, all_consensus

            if len(superfamily_annotations) == 2:
                if all_relation == "coincidencia_exacta":
                    rule = "ConoSorter+2DB_coinciden"
                else:
                    rule = "ConoSorter+2DB_compatibles_" + suffix
                return "ALTA", rule, all_consensus

            if len(superfamily_annotations) == 3:
                if all_relation == "coincidencia_exacta":
                    rule = "ConoSorter+3DB_coinciden"
                else:
                    rule = "ConoSorter+3DB_compatibles_" + suffix
                return "ALTA", rule, all_consensus

        # Regla histórica del proyecto: si ConoSorter y al menos dos
        # bases independientes apoyan el mismo consenso, conservar ALTA
        # aunque una tercera base de similitud sea incompatible.
        #
        # No se aplica a M2/O4 cuando existe una familia externa
        # incompatible, porque esos casos ya fueron detenidos por
        # HARD_PROJECT_CONFLICT_RELATIONS.
        majority_counts = {}
        majority_relations = {}

        for is_compatible, pair_consensus, pair_relation in pair_results:
            if not is_compatible or not pair_consensus:
                continue
            majority_counts[pair_consensus] = (
                majority_counts.get(pair_consensus, 0) + 1
            )
            majority_relations.setdefault(
                pair_consensus,
                set(),
            ).add(pair_relation)

        majority_candidates = [
            (consensus, count)
            for consensus, count in majority_counts.items()
            if count >= 2
        ]

        if len(majority_candidates) == 1:
            majority_consensus, majority_count = majority_candidates[0]
            relations = majority_relations[majority_consensus]

            if relations == {"coincidencia_exacta"}:
                rule = (
                    f"ConoSorter+{majority_count}DB_coinciden_"
                    "1DB_en_conflicto"
                )
            else:
                suffixes = sorted(
                    relation_suffix(relation)
                    for relation in relations
                )
                rule = (
                    f"ConoSorter+{majority_count}DB_compatibles_"
                    + "_".join(suffixes)
                    + "_1DB_en_conflicto"
                )

            return "ALTA", rule, majority_consensus

        # Solo una base es compatible con ConoSorter; las demás difieren.
        if n_compatible == 1:
            _, pair_consensus, pair_relation = compatible_pairs[0]
            suffix = relation_suffix(pair_relation)
            if pair_relation == "coincidencia_exacta":
                rule = "ConoSorter+1DB_coinciden"
            else:
                rule = (
                    "ConoSorter+1DB_compatible_otras_difieren_"
                    + suffix
                )
            return "MEDIA", rule, pair_consensus or conosorter

        # Las bases son compatibles entre sí, pero ConoSorter difiere.
        if db_consensus and len(superfamily_annotations) >= 2:
            suffix = relation_suffix(db_relation)
            n_db = len(superfamily_annotations)

            if db_relation == "coincidencia_exacta":
                rule = f"{n_db}DB_consenso_ConoSorter_difiere"
            else:
                rule = (
                    f"{n_db}DB_compatibles_ConoSorter_difiere_"
                    + suffix
                )

            return "MEDIA", rule, db_consensus

        # Dos etiquetas pueden ser compatibles individualmente con una
        # etiqueta genérica, pero incompatibles entre sí (p. ej., I4 e I5).
        if n_compatible >= 2:
            return (
                "BAJA",
                "Compatibilidad_individual_pero_conflicto_entre_fuentes",
                "",
            )

        return "BAJA", "Todas_fuentes_diferentes", ""

    # Secuencias negativas: solo se evalúan las tres bases de similitud.
    if not superfamily_annotations:
        return "SIN_CLASIFICAR", "Sin_anotaciones", ""

    if db_consensus:
        n_db = len(superfamily_annotations)
        suffix = relation_suffix(db_relation)

        if n_db >= 2:
            if db_relation == "coincidencia_exacta":
                rule = f"{n_db}DB_consenso"
            else:
                rule = f"{n_db}DB_compatibles_{suffix}"
            return "MEDIA", rule, db_consensus

        return "BAJA", "1DB_unica", db_consensus

    return "BAJA", "DB_en_conflicto", ""



def classify_sequence(
    is_positive: bool,
    conosorter: str,
    conoserver: str,
    ncbi: str,
    uniprot: str,
) -> Tuple[str, str, str]:
    """
    Clasifica conservando las etiquetas originales.

    Las clases peptídicas se traducen internamente a su superfamilia génica
    únicamente para evaluar compatibilidad y consenso.
    """
    raw_values = {
        "Signal": conosorter,
        "ConoServer": conoserver,
        "NCBI": ncbi,
        "UniProt": uniprot,
    }

    canonical = {}
    mapping_notes = []

    for source, raw_value in raw_values.items():
        normalized, note = normalize_superfamily_for_compare_with_note(
            raw_value
        )
        canonical[source] = normalized
        if note:
            mapping_notes.append(f"{source}_{note}")

    certidumbre, regla, consenso = _classify_canonical_sequence(
        is_positive=is_positive,
        conosorter=canonical["Signal"],
        conoserver=canonical["ConoServer"],
        ncbi=canonical["NCBI"],
        uniprot=canonical["UniProt"],
    )

    if mapping_notes:
        regla = (
            regla
            + "__mapeo_nomenclatural_"
            + "__".join(mapping_notes)
        )

    return certidumbre, regla, consenso
