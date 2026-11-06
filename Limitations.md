
## 1. Design notes and known limitations

**`filtered_contigs.fasta` is the contract, and it is not always filtered.**
Stage 02 always reads that file. It holds TransRate's `good.*.fa` for
paired-end groups and the **raw** assembly for single-end ones, because
TransRate's read-metrics mode requires both mates. The `filter_method` column
in the assembly manifest records which happened. **Filtered and unfiltered
contig sets are not comparable — stratify on that column** before comparing
anything across samples.

**A missing container runtime fails loudly, on purpose.** Exit codes 10/11/12
mean Docker is stopped or Apptainer is not loaded — a setup fault fixed in
seconds. Any other exit code degrades gracefully and is recorded, because that
means the tool ran but could not score *this* assembly, which is a legitimate
data outcome.

**Trinity on Apple Silicon is slow.** Its image is `linux/amd64` and runs
emulated; a 9 M-pair library can take many hours. Disable it
(`assembly.assemblers.trinity.enabled: false`) to iterate locally and enable it
for production runs on Linux. The choice propagates through the whole DAG,
including the cross-assembler dereplication.

**NCBI contributes few superfamily labels.** Its feature table rarely states a
superfamily explicitly (most records read `product conotoxin ...` without one),
so in practice the medium-confidence consensus rests on ConoServer and UniProt.
This is a property of the source, not a defect, and is reported per group in
the manifest.

**Reads outside the working directory need an explicit bind.** With
`preprocessing.fastp: false`, Stage 01 reads the raw FASTQs directly. If they
live outside the workdir — reads on `/scratch`, run launched from `/home` — add
the parent directory to `containers.extra_binds`.

**Start-codon refinement.** `annotation.predict.no_refine_starts` is `false`
here. Real assemblies contain far more contigs than TransDecoder's start-codon
PWM needs to train, and refinement sharpens exactly the region this workflow is
most sensitive to: the ORF N-terminus *is* the signal peptide fed to SignalP.
Set it to `true` if you run on very sparse assemblies, where the PWM cannot
train and TransDecoder aborts — but set it for the **whole** dataset, never
per-sample, or ORF sets stop being comparable.

---

## 2. FAIR compliance

| Principle | Implementation |
|---|---|
| **Findable** | `CITATION.cff` renders a "Cite this repository" button; archive a release on Zenodo for a persistent DOI; descriptive keywords in the citation metadata |
| **Accessible** | Public repository, MIT licence, no authentication needed; all tool dependencies resolved from public conda channels and container registries |
| **Interoperable** | Standard formats throughout (TSV, FASTA, FASTQ, XLSX); every stage exposes a documented flat-table manifest; Snakemake is a community standard workflow language |
| **Reusable** | MIT licence with third-party terms enumerated; pinned tool versions in `workflow/envs/*.yaml`; SHA-256 sidecars on every manifest; per-sample summaries record every parameter, tool version and container image that could change a result |

**Provenance.** Every stage emits a manifest recording the parameters and tool
versions used. The annotation logic under `workflow/scripts/lib/` is copied
**verbatim** from a published reference implementation — do not edit it in
place; amend the source and the copy together so provenance stays auditable.

---

---

# 🇪🇸 Español

## 1. Notas de diseño y limitaciones conocidas

**`filtered_contigs.fasta` es el contrato, y no siempre está filtrado.**
Stage 02 siempre lee ese archivo. Contiene el `good.*.fa` de TransRate para
grupos paired-end y el ensamble **crudo** para los single-end, porque el modo
read-metrics de TransRate requiere ambas mates. La columna `filter_method` del
manifiesto de ensamble registra cuál ocurrió. **Los conjuntos filtrados y sin
filtrar no son comparables — estratifica por esa columna** antes de comparar
nada entre muestras.

**La ausencia de runtime de contenedores falla ruidosamente, a propósito.**
Los códigos de salida 10/11/12 significan que Docker está detenido o que
Apptainer no está cargado — una falla de configuración que se arregla en
segundos. Cualquier otro código degrada con gracia y queda registrado, porque
eso significa que la herramienta corrió pero no pudo evaluar *este* ensamble,
lo cual es un desenlace legítimo de los datos.

**Trinity en Apple Silicon es lento.** Su imagen es `linux/amd64` y corre
emulada; una librería de 9 M de pares puede tardar muchas horas. Desactívalo
(`assembly.assemblers.trinity.enabled: false`) para iterar localmente y
actívalo para las corridas de producción en Linux. La decisión se propaga por
todo el DAG, incluida la dereplicación entre ensambladores.

**NCBI aporta pocas etiquetas de superfamilia.** Su tabla de características
rara vez declara la superfamilia explícitamente (la mayoría de los registros
dicen `product conotoxin ...` sin ella), así que en la práctica el consenso de
certidumbre media descansa en ConoServer y UniProt. Es una propiedad de la
fuente, no un defecto, y se reporta por grupo en el manifiesto.

**Las lecturas fuera del directorio de trabajo requieren un bind explícito.**
Con `preprocessing.fastp: false`, Stage 01 lee los FASTQ crudos directamente.
Si viven fuera del workdir — lecturas en `/scratch`, corrida lanzada desde
`/home` — agrega el directorio padre a `containers.extra_binds`.

**Refinamiento del codón de inicio.** `annotation.predict.no_refine_starts` es
`false` aquí. Los ensambles reales contienen muchos más contigs de los que
necesita la PWM de codón de inicio de TransDecoder para entrenar, y el
refinamiento afina justamente la región a la que este flujo es más sensible: el
extremo N-terminal del ORF *es* el péptido señal que se le pasa a SignalP.
Ponlo en `true` si corres sobre ensambles muy escasos, donde la PWM no puede
entrenar y TransDecoder aborta — pero ponlo para **todo** el conjunto de datos,
nunca por muestra, o los conjuntos de ORFs dejan de ser comparables.

---

## 2. Cumplimiento FAIR

| Principio | Implementación |
|---|---|
| **Localizable** | `CITATION.cff` genera el botón "Cite this repository"; archiva un release en Zenodo para obtener un DOI persistente; palabras clave descriptivas en los metadatos de citación |
| **Accesible** | Repositorio público, licencia MIT, sin autenticación; todas las dependencias se resuelven desde canales conda y registros de contenedores públicos |
| **Interoperable** | Formatos estándar en todo el flujo (TSV, FASTA, FASTQ, XLSX); cada etapa expone un manifiesto tabular plano documentado; Snakemake es un lenguaje de flujos estándar en la comunidad |
| **Reutilizable** | Licencia MIT con los términos de terceros enumerados; versiones de herramientas fijadas en `workflow/envs/*.yaml`; archivos SHA-256 para cada manifiesto; los resúmenes por muestra registran cada parámetro, versión de herramienta e imagen de contenedor que podría alterar un resultado |

**Procedencia.** Cada etapa emite un manifiesto que registra los parámetros y
versiones de herramientas usadas. La lógica de anotación bajo
`workflow/scripts/lib/` es una copia **literal** de una implementación de
referencia publicada — no la edites en el lugar; modifica la fuente y la copia
juntas para que la procedencia siga siendo auditable.

---

