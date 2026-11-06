# Bioinformatic workflow for conotoxin discovery and annotation

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A reproducible, portable Snakemake workflow that discovers and annotates
conotoxin precursors from real single- or paired-end venom-duct RNA-seq.
Runs unchanged on **macOS with Docker** and on **Linux HPC under SLURM with
Apptainer**.

**🇬🇧 English below · [🇪🇸 Español más abajo](#-español)**

---

## Background

Venom-duct transcriptome assembly is the standard route to the species-specific
sequence databases that proteotranscriptomic venom characterisation depends on.
Compared with proteomics alone it is faster and cheaper, and it targets the
coding sequences actually expressed in the venom-producing tissue. Short-read
RNA-seq dominates this application because it delivers high coverage at modest
cost — but sequencing artefacts, library-design choices and sample complexity
all propagate into transcript reconstruction and, from there, into annotation.

Conotoxin precursors are a particularly hard reconstruction target. Closely
related precursors pair a **conserved signal region** with a **hypervariable
mature region**, so an assembler must resolve near-identical 5′ ends that
diverge sharply downstream. Assemblers differ in how they handle this, and each
one misses precursors the others recover.

This workflow is the applied outcome of a controlled comparison: de novo
assemblers were benchmarked on simulated reads with known ground truth, and a
confidence-based strategy was developed to grade the resulting candidates.
Two design choices follow directly from that work.

**Three assemblers, then dereplication.** Trinity, rnaSPAdes and RNA-Bloom are
run in parallel on the same reads and their predicted peptides are pooled and
dereplicated. The union recovers precursors that any single assembler misses,
while dereplication keeps the redundancy from inflating the candidate count.

**Four evidence sources, then tiers.** Rather than accepting or rejecting a
candidate outright, the workflow records agreement among the ConoSorter
signal-region assignment and normalised annotations from ConoServer, NCBI and
UniProt, and assigns each candidate to a confidence tier. Sequences with no
meaningful agreement are retained as unclassified rather than discarded, so the
evidence behind every call stays inspectable.

---

## Overview

```
00_preprocessing   fastp per library  →  pooling per group
01_assembly        Trinity + rnaSPAdes + RNA-Bloom  →  TransRate filtering
02_annotation      ORF prediction → dereplication → ConoSorter → SignalP 6
                   → DIAMOND vs 3 Conidae databases → label normalisation
                   → confidence tiers (HIGH / MEDIUM / LOW / UNCLASSIFIED)
```

Candidates are graded by the level of agreement among **four independent
evidence sources**: the ConoSorter signal-region superfamily assignment plus
normalised annotations from ConoServer, NCBI and UniProt.

| Tier | Criterion |
|---|---|
| `ALTA` (high) | ConoSorter-assigned **and** ≥1 external database agrees on the same superfamily |
| `MEDIA` (medium) | Partial support, or consensus among external databases despite ConoSorter disagreement |
| `BAJA` (low) | Only weak or conflicting support |
| `SIN_CLASIFICAR` (unclassified) | No meaningful agreement |

---

## 1. The two units of analysis

This distinction shapes the entire layout and is the most important thing to
understand before writing your sample sheet.

| | Unit | Rationale |
|---|---|---|
| Stage 00 | **library** (`sample_id`) | Adapter content and base quality differ between lanes and sequencing runs. Each library is trimmed with its own adapter detection and keeps an individual QC report — which is also how a failed replicate becomes visible instead of being averaged away. |
| Stages 01–02 | **group** | Replicates (or split files) assembled separately fragment transcripts that the pooled depth recovers intact. This matters disproportionately for conotoxins: precursors are short (~60–150 aa) and expression is highly uneven, so a partial ORF from one replicate typically fails the signal-peptide test downstream and disappears from the annotation entirely. |

A library with an empty `group` forms a group of one, reproducing
one-assembly-per-library behaviour.

---

## 2. Requirements

### Software

| Requirement | Version | Notes |
|---|---|---|
| Snakemake | ≥ 9.0 | `pip install snakemake` |
| Conda / Mamba | any recent | Environments are created automatically |
| Docker **or** Apptainer/Singularity | any | For Trinity, TransRate and the SignalP fallback |
| SLURM executor plugin | — | HPC only: `pip install snakemake-executor-plugin-slurm` |

All other tools (fastp, rnaSPAdes, RNA-Bloom, TransDecoder, DIAMOND, MMseqs2,
HMMER, gawk) are installed automatically into conda environments.

### Data you must provide

| Resource | Config key | How to obtain |
|---|---|---|
| Conidae protein FASTAs + DIAMOND indexes for **ConoServer**, **NCBI** and **UniProt** | `annotation.databases.root` | Download from each provider; build indexes with `diamond makedb` |
| NCBI feature table (`.txt`) and UniProt labels (`.xlsx`) | same | Exported alongside the FASTAs |
| ConoSorter v1.1 ZIP | `annotation.conosorter.local_zip` | [SourceForge](https://sourceforge.net/projects/conosorter/) (GPLv3) |
| SignalP 6.0 | `annotation.signalp6.*` | See below |

> **SignalP 6.0 is not redistributable.** Obtain it under the free academic
> licence from [DTU Health Tech](https://services.healthtech.dtu.dk/services/SignalP-6.0/)
> and place the tarball at the configured path. If no tarball is found the
> workflow automatically falls back to a public container image; the route
> actually used is recorded as `signalp_route` in the output manifest, because
> the container is third-party and not version-pinned by this project.

### Hardware

| Stage | Typical requirement |
|---|---|
| fastp | 8 GB RAM |
| Trinity | 32–64 GB RAM, many hours (native Linux strongly preferred) |
| rnaSPAdes | 32–64 GB RAM |
| TransRate | 16 GB RAM **+ ~7 GB free disk per 9 M read pairs** (it cannot read gzip, so both mates are decompressed for the duration of the rule) |
| SignalP 6 | 12–24 GB RAM (loads a ~1.5 GB transformer) |

---

## 3. Installation

```bash
git clone https://github.com/rjegr/conotoxin-discovery-workflow.git
cd conotoxin-discovery-workflow
```

### 3.1 Create `resources/`

The repository ships **code only**. Databases and licensed software are not
redistributable, so you must supply them under a `resources/` directory at the
repository root. It is gitignored, so nothing you put there will ever be
committed by accident.

```
resources/
├── conidae_dbs/20260701/
│   ├── conoserver/
│   │   ├── conoserver_240620_protein.fa
│   │   └── conoserver.dmnd
│   ├── NCBI/proteinas/
│   │   ├── NCBI_Conidae_proteins_26_6_26.fasta
│   │   └── NCBI_Conidae_proteins_26_6_26.dmnd
│   ├── UniProt/
│   │   ├── Uniprot_Conidae_proteins_26_6_26.fasta
│   │   └── Uniprot_Conidae_proteins_26_6_26.dmnd
│   ├── base_datos_anotacion_ncbi.txt        # NCBI feature table
│   └── base_datos_anotacion_uniprot.xlsx    # UniProt "Protein names" table
└── software/
    ├── conosorter/ConoSorter_v1.1.zip
    └── signalp6/signalp-6.0h.fast.tar.gz    # optional; container used if absent
```

Three ways to provide it, all equivalent:

```bash
# a) copy a prepared snapshot in
cp -R /path/to/prepared/resources resources

# b) symlink a shared copy (keeps one copy per machine)
ln -s /shared/storage/conotoxin-resources resources

# c) point the config at an absolute path instead — usual on HPC, where the
#    databases live on shared storage
#    config/config.yaml:
#      annotation.databases.root: "/shared/db/conidae_dbs/20260701"
```

Build the DIAMOND indexes if you only have the FASTAs:

```bash
for db in conoserver/conoserver_240620_protein \
          NCBI/proteinas/NCBI_Conidae_proteins_26_6_26 \
          UniProt/Uniprot_Conidae_proteins_26_6_26 ; do
    diamond makedb --in resources/conidae_dbs/20260701/${db}.fa* \
                   -d resources/conidae_dbs/20260701/${db}
done
```

The exact filenames are configurable — see `annotation.databases` in
`config/config.yaml` — so an existing snapshot with different names needs a
config edit, not a rename.

### 3.2 Verify

```bash
# Static validation — no data required. Checks syntax, YAMLs, vendored
# scripts, the sample sheet and the container runtime.
bash scripts_check.sh

# Resource validation — confirms every database, index and licensed tool is
# where the config says. Run this BEFORE spending CPU on assembly.
snakemake --profile profiles/local results/run01/qc/annotation_preflight.txt
cat results/run01/qc/annotation_preflight.txt
```

The preflight also reports which SignalP route it will take (`tarball` or
`container`), so you know before the run starts.

---

## 4. Input: `config/samples.tsv`

Tab-separated. `group` and `fastq_2` are optional columns.

```tsv
sample_id	group	fastq_1	fastq_2
Cmar_rep1	Cmarmoreus_VD	/abs/path/Cmar_rep1_R1.fastq.gz	/abs/path/Cmar_rep1_R2.fastq.gz
Cmar_rep2	Cmarmoreus_VD	/abs/path/Cmar_rep2_R1.fastq.gz	/abs/path/Cmar_rep2_R2.fastq.gz
Cpurp_SE	Cpurpurascens	/abs/path/Cpurp_SE.fastq.gz
```

Rules enforced by the validator:

- Paths must be **absolute** — the workflow may resolve them from a different
  directory on a compute node than on the submit node.
- `sample_id` and `group` may contain only `[A-Za-z0-9._-]`; they become
  directory names and Snakemake wildcards.
- An empty `fastq_2` marks a single-end library.
- **All libraries in a group must share one layout.** Concatenating unpaired
  with paired reads produces an incoherent read set that assemblers either
  reject or silently mis-assemble.

> ⚠️ Write it with **real tab characters**. A sheet indented with spaces parses
> as a single column. Verify with `cat -A config/samples.tsv` — tabs appear
> as `^I`.

---

## 5. Configuration

Edit `config/config.yaml`. The keys you are most likely to change:

```yaml
run_id: "run01"                       # results/{run_id}/...
samples: "config/samples.tsv"

preprocessing:
  fastp: true                         # false → assemble raw reads

assembly:
  assemblers:
    trinity:   {enabled: true,  max_memory_gb: 32, threads: 8}
    rnaspades: {enabled: true,  max_memory_gb: 32, threads: 8}
    rnabloom:  {enabled: true,  max_memory_gb: 16, threads: 8}
  transrate:
    enabled: true

annotation:
  databases:
    root: "resources/conidae_dbs/20260701"
  predict:
    no_refine_starts: false           # see note below
```

> **Docker Desktop's memory is not your machine's memory.** Check with
> `docker info | grep -i "total memory"` and keep `max_memory_gb` at roughly
> 75% of that value. Setting it at or above the container limit invites a
> silent OOM kill.

---

## 6. Running

```bash
# 1. Validate the sample sheet — the only step allowed to hard-fail on input
snakemake --profile profiles/local validate
cat results/run01/qc/samples_validation_report.txt
cat results/run01/qc/groups_validated.tsv

# 2. Verify annotation resources BEFORE spending CPU on assembly
snakemake --profile profiles/local results/run01/qc/annotation_preflight.txt
cat results/run01/qc/annotation_preflight.txt

# 3. Run — stage by stage, or all at once
snakemake --profile profiles/local preprocessing
snakemake --profile profiles/local assembly
snakemake --profile profiles/local annotation
snakemake --profile profiles/local all
```

On HPC, swap the profile. **Nothing else changes** — the container helper
detects Apptainer and translates the flags automatically:

```bash
snakemake --profile profiles/slurm all
```

Before first use on a cluster, edit `profiles/slurm/config.yaml`:
`slurm_partition`, `slurm_account` (if required), and `conda-prefix` (point it
at shared storage, **not** `$HOME` if you have a small quota).

> Always dry-run first (`-n`) and read the **Job stats** table. If it lists
> rules you did not expect, stop and investigate before committing cluster
> time.

---

## 7. Outputs

```
results/{run_id}/
├── 00_preprocessing/
│   ├── {sample}/                     clean reads, fastp.html/.json, summary
│   ├── _groups/{group}/              pooled reads + concat summary
│   └── preprocessing_manifest.tsv    one row per LIBRARY
├── 01_assembly/
│   ├── {group}/{asm}/                contigs.fasta, filtered_contigs.fasta
│   └── assembly_manifest.tsv         one row per (group, assembler)
├── 02_annotation/
│   ├── {group}/{asm}/predict/        TransDecoder ORFs
│   ├── {group}/dedup/                dereplicated peptides
│   ├── {group}/conosorter/           RegEx + pHMM candidates
│   ├── {group}/signalp/              SP-positive, split assigned/unassigned
│   ├── {group}/annotate/             DIAMOND hits + resolved labels
│   ├── {group}/classify/
│   │   └── classification.xlsx       ← THE BIOLOGICAL DELIVERABLE
│   └── annotation_manifest.tsv       one row per GROUP
└── qc/
    ├── samples_validation_report.txt
    ├── annotation_preflight.txt
    ├── annotation_invariants.txt
    └── pipeline_gate.OK              written only on a clean run
```

`classification.xlsx` contains five sheets: one per confidence tier plus all
sequences.

Every manifest ships with a `.sha256` sidecar.

---

## 8. Reference results

End-to-end run on a *Conus* venom-duct library (3 pooled paired-end libraries,
9,324,984 read pairs, 2 assemblers):

| Step | Count |
|---|---|
| ORFs predicted | 14,227 |
| After dereplication | 9,075 (64%) |
| ConoSorter candidates | 8,836 (97%) |
| **SignalP-positive** | **597 (6.8%)** |
| ConoSorter-assigned / unassigned | 213 / 384 |
| **HIGH confidence** | **170** |
| MEDIUM / LOW / UNCLASSIFIED | 49 / 93 / 285 |

**Interpreting the SignalP rate.** 6.8% is the expected order of magnitude for
a real transcriptome: a venom duct expresses mostly non-toxin transcripts, and
the secretion-signal screen is what removes them. (On simulated data derived
from curated conotoxins the same step retains ~75%, because almost every input
is genuinely a toxin.) 170 high-confidence candidates is biologically
plausible — *Conus* species typically express 100–400 conotoxins.

---


## 9. Tool Reference

Please cite this workflow (see `CITATION.cff`) **and** the tools it
orchestrates:

- **fastp** — Chen S. *et al.* (2018) *Bioinformatics* 34:i884–i890
- **Trinity** — Grabherr M.G. *et al.* (2011) *Nat Biotechnol* 29:644–652
- **rnaSPAdes** — Bushmanova E. *et al.* (2019) *GigaScience* 8:giz100
- **RNA-Bloom** — Nip K.M. *et al.* (2020) *Genome Res* 30:1191–1200
- **TransRate** — Smith-Unna R. *et al.* (2016) *Genome Res* 26:1134–1144
- **TransDecoder** — Haas B.J. *et al.* (2013) *Nat Protoc* 8:1494–1512
- **DIAMOND** — Buchfink B. *et al.* (2021) *Nat Methods* 18:366–368
- **MMseqs2** — Steinegger M. & Söding J. (2017) *Nat Biotechnol* 35:1026–1028
- **ConoSorter** — Lavergne V. *et al.* (2013) *BMC Genomics* 14:708
- **SignalP 6.0** — Teufel F. *et al.* (2022) *Nat Biotechnol* 40:1023–1025
- **ConoServer** — Kaas Q. *et al.* (2012) *Bioinformatics* 28:1461–1462
- **Snakemake** — Mölder F. *et al.* (2021) *F1000Research* 10:33

---


## 10. Licence

MIT — see [LICENSE](LICENSE). The MIT terms cover **only the workflow code in
this repository**. Orchestrated tools and reference databases carry their own
licences, enumerated in `LICENSE`. SignalP 6.0 in particular is **not
redistributable** and must be obtained by each user under the DTU academic
licence.

---
---

# 🇪🇸 Español

# Flujo de trabajo bioinformático para el descubrimiento y anotación de conotoxinas

Flujo de trabajo reproducible y portable en Snakemake que descubre y anota
precursores de conotoxinas a partir de RNA-seq real de ducto venenoso, en
formato single-end o paired-end. Corre sin modificaciones en **macOS con
Docker** y en **HPC Linux bajo SLURM con Apptainer**.

## Antecedentes

El ensamble del transcriptoma de ducto venenoso es la vía estándar para
construir las bases de datos específicas de especie de las que depende la
caracterización proteotranscriptómica de venenos. Frente a la proteómica sola
es más rápido y económico, y apunta directamente a las secuencias codificantes
que realmente se expresan en el tejido productor de veneno. El RNA-seq de
lecturas cortas domina esta aplicación porque ofrece alta cobertura a costo
moderado — pero los artefactos de secuenciación, las decisiones de diseño de
librería y la complejidad de la muestra se propagan hacia la reconstrucción de
transcritos y, desde ahí, hacia la anotación.

Los precursores de conotoxinas son un blanco de reconstrucción especialmente
difícil. Precursores estrechamente relacionados combinan una **región señal
conservada** con una **región madura hipervariable**, de modo que el ensamblador
debe resolver extremos 5′ casi idénticos que divergen abruptamente río abajo.
Los ensambladores difieren en cómo manejan esto, y cada uno pierde precursores
que los otros sí recuperan.

Este flujo de trabajo es el resultado aplicado de una comparación controlada:
se evaluaron ensambladores de novo sobre lecturas simuladas con verdad de
fondo conocida, y se desarrolló una estrategia basada en certidumbre para
graduar los candidatos resultantes. De ese trabajo se desprenden dos decisiones
de diseño.

**Tres ensambladores, luego dereplicación.** Trinity, rnaSPAdes y RNA-Bloom
corren en paralelo sobre las mismas lecturas y sus péptidos predichos se
combinan y dereplican. La unión recupera precursores que cualquier ensamblador
individual pierde, mientras que la dereplicación evita que la redundancia infle
el conteo de candidatos.

**Cuatro fuentes de evidencia, luego niveles.** En lugar de aceptar o rechazar
un candidato de tajo, el flujo registra el acuerdo entre la asignación de
ConoSorter en la región señal y las anotaciones normalizadas de ConoServer,
NCBI y UniProt, y asigna cada candidato a un nivel de certidumbre. Las
secuencias sin acuerdo significativo se conservan como no clasificadas en lugar
de descartarse, de modo que la evidencia detrás de cada decisión permanece
inspeccionable.

## Panorama general

```
00_preprocessing   fastp por librería  →  agrupamiento por grupo
01_assembly        Trinity + rnaSPAdes + RNA-Bloom  →  filtrado TransRate
02_annotation      predicción de ORFs → dereplicación → ConoSorter → SignalP 6
                   → DIAMOND vs 3 bases Conidae → normalización de etiquetas
                   → niveles de certidumbre (ALTA / MEDIA / BAJA / SIN_CLASIFICAR)
```

Los candidatos se gradúan por el nivel de acuerdo entre **cuatro fuentes de
evidencia independientes**: la asignación de superfamilia de ConoSorter (región
señal) más las anotaciones normalizadas de ConoServer, NCBI y UniProt.

| Nivel | Criterio |
|---|---|
| `ALTA` | Asignado por ConoSorter **y** ≥1 base externa coincide en la misma superfamilia |
| `MEDIA` | Apoyo parcial, o consenso entre bases externas pese a discrepancia con ConoSorter |
| `BAJA` | Solo apoyo débil o contradictorio |
| `SIN_CLASIFICAR` | Sin acuerdo significativo |

---

## 1. Las dos unidades de análisis

Esta distinción define toda la arquitectura y es lo más importante que hay que
entender antes de escribir el manifiesto de muestras.

| | Unidad | Razón |
|---|---|---|
| Stage 00 | **librería** (`sample_id`) | El contenido de adaptadores y la calidad de bases difieren entre lanes y corridas de secuenciación. Cada librería se limpia con su propia detección de adaptadores y conserva un reporte de QC individual — que además es la forma de detectar una réplica fallida en lugar de promediarla. |
| Stages 01–02 | **grupo** | Ensamblar réplicas (o archivos partidos) por separado fragmenta transcritos que la profundidad combinada recupera completos. Esto pesa desproporcionadamente en conotoxinas: los precursores son cortos (~60–150 aa) y su expresión es muy desigual, así que un ORF parcial de una réplica típicamente falla la prueba de péptido señal río abajo y desaparece por completo de la anotación. |

Una librería con `group` vacío forma un grupo de una sola, reproduciendo el
comportamiento de un ensamble por librería.

---

## 2. Requisitos

### Software

| Requisito | Versión | Notas |
|---|---|---|
| Snakemake | ≥ 9.0 | `pip install snakemake` |
| Conda / Mamba | reciente | Los ambientes se crean automáticamente |
| Docker **o** Apptainer/Singularity | cualquiera | Para Trinity, TransRate y el respaldo de SignalP |
| Plugin ejecutor SLURM | — | Solo HPC: `pip install snakemake-executor-plugin-slurm` |

Las demás herramientas (fastp, rnaSPAdes, RNA-Bloom, TransDecoder, DIAMOND,
MMseqs2, HMMER, gawk) se instalan solas en ambientes conda.

### Datos que debes proporcionar

| Recurso | Clave de config | Cómo obtenerlo |
|---|---|---|
| FASTAs de proteínas Conidae + índices DIAMOND de **ConoServer**, **NCBI** y **UniProt** | `annotation.databases.root` | Descargar de cada proveedor; construir índices con `diamond makedb` |
| Tabla de características NCBI (`.txt`) y etiquetas UniProt (`.xlsx`) | igual | Se exportan junto con los FASTAs |
| ZIP de ConoSorter v1.1 | `annotation.conosorter.local_zip` | [SourceForge](https://sourceforge.net/projects/conosorter/) (GPLv3) |
| SignalP 6.0 | `annotation.signalp6.*` | Ver abajo |

> **SignalP 6.0 no es redistribuible.** Obtenlo bajo la licencia académica
> gratuita de [DTU Health Tech](https://services.healthtech.dtu.dk/services/SignalP-6.0/)
> y coloca el tarball en la ruta configurada. Si no se encuentra ningún
> tarball, el flujo recurre automáticamente a una imagen de contenedor pública;
> la ruta realmente usada se registra como `signalp_route` en el manifiesto de
> salida, porque el contenedor es de terceros y este proyecto no fija su
> versión.

### Hardware

| Etapa | Requisito típico |
|---|---|
| fastp | 8 GB RAM |
| Trinity | 32–64 GB RAM, varias horas (se recomienda fuertemente Linux nativo) |
| rnaSPAdes | 32–64 GB RAM |
| TransRate | 16 GB RAM **+ ~7 GB de disco libre por cada 9 M de pares** (no lee gzip, así que ambas mates se descomprimen mientras dura la regla) |
| SignalP 6 | 12–24 GB RAM (carga un transformer de ~1.5 GB) |

---

## 3. Instalación

```bash
git clone https://github.com/rjegr/conotoxin-discovery-workflow.git
cd conotoxin-discovery-workflow
```

### 3.1 Crear `resources/`

El repositorio contiene **solo código**. Las bases de datos y el software con
licencia no son redistribuibles, así que debes proporcionarlos en un directorio
`resources/` en la raíz del repositorio. Está en `.gitignore`, así que nada de
lo que pongas ahí se subirá por accidente.

```
resources/
├── conidae_dbs/20260701/
│   ├── conoserver/
│   │   ├── conoserver_240620_protein.fa
│   │   └── conoserver.dmnd
│   ├── NCBI/proteinas/
│   │   ├── NCBI_Conidae_proteins_26_6_26.fasta
│   │   └── NCBI_Conidae_proteins_26_6_26.dmnd
│   ├── UniProt/
│   │   ├── Uniprot_Conidae_proteins_26_6_26.fasta
│   │   └── Uniprot_Conidae_proteins_26_6_26.dmnd
│   ├── base_datos_anotacion_ncbi.txt        # tabla de características NCBI
│   └── base_datos_anotacion_uniprot.xlsx    # tabla "Protein names" de UniProt
└── software/
    ├── conosorter/ConoSorter_v1.1.zip
    └── signalp6/signalp-6.0h.fast.tar.gz    # opcional; si falta se usa contenedor
```

Tres formas de proporcionarlo, todas equivalentes:

```bash
# a) copiar un snapshot ya preparado
cp -R /ruta/al/snapshot/resources resources

# b) symlink a una copia compartida (mantiene una sola copia por máquina)
ln -s /almacenamiento/compartido/conotoxin-resources resources

# c) apuntar el config a una ruta absoluta — lo habitual en HPC, donde las
#    bases viven en almacenamiento compartido
#    config/config.yaml:
#      annotation.databases.root: "/shared/db/conidae_dbs/20260701"
```

Construye los índices DIAMOND si solo tienes los FASTAs:

```bash
for db in conoserver/conoserver_240620_protein \
          NCBI/proteinas/NCBI_Conidae_proteins_26_6_26 \
          UniProt/Uniprot_Conidae_proteins_26_6_26 ; do
    diamond makedb --in resources/conidae_dbs/20260701/${db}.fa* \
                   -d resources/conidae_dbs/20260701/${db}
done
```

Los nombres exactos de archivo son configurables — ver `annotation.databases`
en `config/config.yaml` — así que un snapshot existente con otros nombres
requiere editar el config, no renombrar archivos.

### 3.2 Verificar

```bash
# Validación estática — no requiere datos. Revisa sintaxis, YAMLs, scripts
# vendorizados, el manifiesto de muestras y el runtime de contenedores.
bash scripts_check.sh

# Validación de recursos — confirma que cada base, índice y herramienta con
# licencia está donde el config dice. Córrela ANTES de gastar CPU en ensamblar.
snakemake --profile profiles/local results/run01/qc/annotation_preflight.txt
cat results/run01/qc/annotation_preflight.txt
```

El preflight además reporta qué ruta de SignalP va a tomar (`tarball` o
`container`), así lo sabes antes de que arranque la corrida.

---

## 4. Entrada: `config/samples.tsv`

Separado por tabuladores. `group` y `fastq_2` son columnas opcionales.

```tsv
sample_id	group	fastq_1	fastq_2
Cmar_rep1	Cmarmoreus_VD	/ruta/abs/Cmar_rep1_R1.fastq.gz	/ruta/abs/Cmar_rep1_R2.fastq.gz
Cmar_rep2	Cmarmoreus_VD	/ruta/abs/Cmar_rep2_R1.fastq.gz	/ruta/abs/Cmar_rep2_R2.fastq.gz
Cpurp_SE	Cpurpurascens	/ruta/abs/Cpurp_SE.fastq.gz
```

Reglas que aplica el validador:

- Las rutas deben ser **absolutas** — el flujo puede resolverlas desde un
  directorio distinto en un nodo de cómputo que en el nodo de envío.
- `sample_id` y `group` solo admiten `[A-Za-z0-9._-]`; se convierten en nombres
  de directorio y en comodines de Snakemake.
- Un `fastq_2` vacío indica librería single-end.
- **Todas las librerías de un grupo deben compartir el mismo layout.**
  Concatenar lecturas no pareadas con pareadas produce un conjunto incoherente
  que los ensambladores rechazan o, peor, ensamblan mal en silencio.

> ⚠️ Escríbelo con **tabuladores reales**. Un manifiesto indentado con espacios
> se interpreta como una sola columna. Verifica con
> `cat -A config/samples.tsv` — los tabs aparecen como `^I`.

---

## 5. Configuración

Edita `config/config.yaml`. Las claves que más probablemente cambiarás:

```yaml
run_id: "run01"                       # results/{run_id}/...
samples: "config/samples.tsv"

preprocessing:
  fastp: true                         # false → ensamblar lecturas crudas

assembly:
  assemblers:
    trinity:   {enabled: true,  max_memory_gb: 32, threads: 8}
    rnaspades: {enabled: true,  max_memory_gb: 32, threads: 8}
    rnabloom:  {enabled: true,  max_memory_gb: 16, threads: 8}
  transrate:
    enabled: true

annotation:
  databases:
    root: "resources/conidae_dbs/20260701"
  predict:
    no_refine_starts: false           # ver nota abajo
```

> **La memoria de Docker Desktop no es la de tu máquina.** Revísala con
> `docker info | grep -i "total memory"` y mantén `max_memory_gb` en ~75% de
> ese valor. Ponerlo igual o por encima del límite del contenedor invita a un
> OOM kill silencioso.

---

## 6. Ejecución

```bash
# 1. Validar el manifiesto — el único paso autorizado a fallar duro por entrada
snakemake --profile profiles/local validate
cat results/run01/qc/samples_validation_report.txt
cat results/run01/qc/groups_validated.tsv

# 2. Verificar los recursos de anotación ANTES de gastar CPU en ensamblar
snakemake --profile profiles/local results/run01/qc/annotation_preflight.txt
cat results/run01/qc/annotation_preflight.txt

# 3. Correr — etapa por etapa, o todo de una vez
snakemake --profile profiles/local preprocessing
snakemake --profile profiles/local assembly
snakemake --profile profiles/local annotation
snakemake --profile profiles/local all
```

En HPC, cambia el profile. **Nada más cambia** — el helper de contenedores
detecta Apptainer y traduce las banderas automáticamente:

```bash
snakemake --profile profiles/slurm all
```

Antes del primer uso en un cluster, edita `profiles/slurm/config.yaml`:
`slurm_partition`, `slurm_account` (si tu sitio lo exige) y `conda-prefix`
(apúntalo a almacenamiento compartido, **no** a `$HOME` si tienes cuota chica).

> Siempre haz primero un dry-run (`-n`) y lee la tabla **Job stats**. Si lista
> reglas que no esperabas, detente e investiga antes de comprometer tiempo de
> cluster.

---

## 7. Salidas

```
results/{run_id}/
├── 00_preprocessing/
│   ├── {sample}/                     lecturas limpias, fastp.html/.json, resumen
│   ├── _groups/{group}/              lecturas agrupadas + resumen de concat
│   └── preprocessing_manifest.tsv    una fila por LIBRERÍA
├── 01_assembly/
│   ├── {group}/{asm}/                contigs.fasta, filtered_contigs.fasta
│   └── assembly_manifest.tsv         una fila por (grupo, ensamblador)
├── 02_annotation/
│   ├── {group}/{asm}/predict/        ORFs de TransDecoder
│   ├── {group}/dedup/                péptidos dereplicados
│   ├── {group}/conosorter/           candidatos RegEx + pHMM
│   ├── {group}/signalp/              SP-positivos, divididos assigned/unassigned
│   ├── {group}/annotate/             hits DIAMOND + etiquetas resueltas
│   ├── {group}/classify/
│   │   └── classification.xlsx       ← EL ENTREGABLE BIOLÓGICO
│   └── annotation_manifest.tsv       una fila por GRUPO
└── qc/
    ├── samples_validation_report.txt
    ├── annotation_preflight.txt
    ├── annotation_invariants.txt
    └── pipeline_gate.OK              se escribe solo si la corrida cierra limpia
```

`classification.xlsx` contiene cinco hojas: una por nivel de certidumbre más
todas las secuencias.

Cada manifiesto trae su archivo `.sha256` acompañante.

---

## 8. Resultados de referencia

Corrida end-to-end sobre una librería de ducto venenoso de *Conus* (3 librerías
paired-end agrupadas, 9,324,984 pares de lecturas, 2 ensambladores):

| Paso | Conteo |
|---|---|
| ORFs predichos | 14,227 |
| Tras dereplicación | 9,075 (64%) |
| Candidatos ConoSorter | 8,836 (97%) |
| **SignalP-positivos** | **597 (6.8%)** |
| Asignados / no asignados por ConoSorter | 213 / 384 |
| **Certidumbre ALTA** | **170** |
| MEDIA / BAJA / SIN_CLASIFICAR | 49 / 93 / 285 |

**Cómo interpretar la tasa de SignalP.** 6.8% es el orden de magnitud esperado
para un transcriptoma real: un ducto venenoso expresa mayoritariamente
transcritos que no son toxinas, y el tamizaje por señal de secreción es lo que
los elimina. (Sobre datos simulados derivados de conotoxinas curadas, ese mismo
paso retiene ~75%, porque casi toda la entrada es genuinamente una toxina.) 170
candidatos de alta certidumbre es biológicamente plausible — las especies de
*Conus* típicamente expresan entre 100 y 400 conotoxinas.

---

## 9. Cómo citar

Cita este flujo de trabajo (ver `CITATION.cff`) **y** las herramientas que
orquesta; la lista completa de referencias está en la
[sección en inglés](#11-citing).

---


## 10. Licencia

MIT — ver [LICENSE](LICENSE). Los términos MIT cubren **únicamente el código
del flujo de trabajo en este repositorio**. Las herramientas orquestadas y las
bases de datos de referencia tienen sus propias licencias, enumeradas en
`LICENSE`. SignalP 6.0 en particular **no es redistribuible** y cada usuario
debe obtenerlo bajo la licencia académica de DTU.

---

## Contact:

Ricardo Gómez-Reyes · rgomez41@uabc.edu.mx · 
