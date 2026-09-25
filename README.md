# PhenoGeneR

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![R](https://img.shields.io/badge/R-%3E%3D4.0.0-brightgreen.svg)
![Python](https://img.shields.io/badge/Python-%3E%3D3.8-brightgreen.svg)

**PhenoGeneR** is an open-source R/Python workflow for **source-aware retrieval, gene-symbol harmonisation, evidence comparison and transparent prioritisation** of human genes associated with a **phenotype, disease or trait term**.

The workflow is designed as an upstream evidence-integration layer. It retrieves source-specific gene evidence, preserves provenance, harmonises gene symbols against the NCBI human gene reference, quantifies cross-source heterogeneity and provides interpretable ranking baselines. PhenoGeneR is **not** intended to assign causal genes automatically.

---

## Overview

Phenotype-, disease- and trait-associated gene evidence is distributed across clinical, genetic-association, literature, expression, pathway, protein and interaction resources. These sources differ in terminology, interfaces, evidence semantics and output structure.

PhenoGeneR provides a single workflow to:

1. query multiple biological resources from a phenotype, disease or trait term;
2. retain source-specific evidence and provenance;
3. combine source-level gene sets;
4. validate human gene symbols against the NCBI human gene reference;
5. resolve recognised aliases and historical symbols to current official symbols;
6. retain unresolved identifiers separately for audit;
7. quantify source coverage and cross-source heterogeneity;
8. construct transparent full-evidence and direct-evidence rankings; and
9. reproduce the benchmark analyses reported in the accompanying Bioinformatics Application Note.

The final manuscript benchmark evaluates **12 source modules** across **13 benchmark terms**.

---

## Evaluated data sources

| Source | Script | Evidence/context represented | Direct-evidence ranking |
|---|---|---|---:|
| ClinVar | `clinvar.R` | Curated clinical variant/gene evidence | Yes |
| Human Phenotype Ontology (HPO) | `hpo.R` | Curated phenotype-gene evidence | Yes |
| OMIM | `omim.R` | Mendelian disease-gene evidence | Yes* |
| GWAS Catalog | `gwasrapidd.R` | Human genetic association evidence and mapped genes | Yes |
| Open Targets | `opentargets.R` | Integrated disease-target evidence | Yes |
| PubMed + PubTator3 | `pubmed.R` | Literature-derived gene evidence | No |
| GTEx | `gtex.R` | Expression/eQTL context | No |
| KEGG | `kegg.R` | Pathway context | No |
| Reactome | `reactome_pathways.R` | Pathway context | No |
| STRING-DB | `string_db.R` | Protein-interaction/network expansion | No |
| UniProt | `uniprot.R` | Protein annotation | No |
| Gene Ontology | `gene_ontology.R` | Functional annotation | No |

\*Only official OMIM associations are treated as direct evidence in the ranking analysis.

### Optional DisGeNET module

`disgenet.R` is retained as an **optional/development module**, but DisGeNET is **not part of the final evaluated 12-source benchmark** reported in the manuscript because reproducible programmatic access was not available for the locked benchmark.

If the optional module is used, provide credentials through the environment:

```bash
export DISGENET_API_KEY="YOUR_KEY"
Rscript disgenet.R migraine
```

Do **not** commit API keys or other credentials to this repository.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/MuhammadMuneeb007/PhenoGeneR.git
cd PhenoGeneR
```

### 2. Install R dependencies

```bash
Rscript requirements.R
```

### 3. Install Python dependencies

Using `pip`:

```bash
pip install -r requirements.txt
```

or Conda:

```bash
conda env create -f environment.yml
conda activate gene-analysis
```

### Requirements

- R >= 4.0
- Python >= 3.8
- internet access for source APIs/services
- source-specific access requirements where applicable

---

## Quick start

### Run the coordinated retrieval workflow

```bash
Rscript download_genes.R migraine
```

For a multi-word term:

```bash
Rscript download_genes.R "high cholesterol"
```

To regenerate outputs instead of using existing files:

```bash
Rscript download_genes.R migraine --force
```

### Run individual source modules

```bash
Rscript clinvar.R migraine
Rscript hpo.R migraine
Rscript omim.R migraine
Rscript gwasrapidd.R migraine
Rscript opentargets.R migraine
Rscript pubmed.R migraine
Rscript gtex.R migraine
Rscript kegg.R migraine
Rscript reactome_pathways.R migraine
Rscript string_db.R migraine
Rscript uniprot.R migraine
Rscript gene_ontology.R migraine
```

---

## Source-specific retrieval

### ClinVar

```bash
Rscript clinvar.R <term>
```

Searches ClinVar and parses structured VCV XML fields to retrieve gene symbols and associated clinical metadata.

### Human Phenotype Ontology

```bash
Rscript hpo.R <term>
```

Uses HPO ontology/annotation resources to identify phenotype-related records and associated genes.

### OMIM

```bash
Rscript omim.R <term>
```

Retrieves Mendelian disease/gene evidence using available OMIM access routes. Only official OMIM associations are treated as direct evidence in the manuscript ranking analysis.

### GWAS Catalog

```bash
Rscript gwasrapidd.R <term>
```

Queries EFO and reported-trait information and maps retrieved variant contexts to genes.

### Open Targets

```bash
Rscript opentargets.R <term>
```

Uses GraphQL disease search and associated-target retrieval.

### PubMed + PubTator3

```bash
Rscript pubmed.R <term>
```

Retrieves term-linked PubMed records and extracts gene annotations using PubTator3.

### GTEx

```bash
Rscript gtex.R <term>
```

Provides expression/eQTL context using phenotype-informed tissue prioritisation and GTEx eGene information.

### KEGG

```bash
Rscript kegg.R <term>
```

Scores human pathway titles against the input term, retrieves matched pathways using `KEGGREST::keggGet()` and parses the pathway `GENE` field as alternating gene-ID/description pairs to obtain gene symbols directly.

### Reactome

```bash
Rscript reactome_pathways.R <term>
```

Matches the input term to Reactome pathway names and maps matched pathways to human gene symbols.

### STRING-DB

```bash
Rscript string_db.R <term>
```

Uses term-resolved seed proteins and expands them through STRING interaction partners. STRING is treated as network-context evidence rather than direct phenotype-gene evidence.

### UniProt

```bash
Rscript uniprot.R <term>
```

Retrieves reviewed human protein records matching query templates and extracts associated gene names and synonyms.

### Gene Ontology

```bash
Rscript gene_ontology.R <term>
```

Retrieves functional annotation context and is treated as indirect functional evidence in downstream ranking analyses.

---

## Output structure

Source-level files are written to:

```text
AllPackagesGenes/
```

Typical outputs include:

```text
migraine_clinvar.csv
migraine_clinvar_genes.csv
migraine_hpo.csv
migraine_hpo_genes.csv
migraine_omim.csv
migraine_omim_genes.csv
migraine_gwasrapidd.csv
migraine_gwasrapidd_genes.csv
migraine_opentargets.csv
migraine_opentargets_genes.csv
migraine_pubmed_pubtator.csv
migraine_pubmed_genes.csv
migraine_gtex_genes.csv
migraine_kegg.csv
migraine_kegg_genes.csv
migraine_reactome_pathways.csv
migraine_reactome_pathways_genes.csv
migraine_string_db.csv
migraine_string_db_genes.csv
migraine_uniprot.csv
migraine_uniprot_genes.csv
migraine_gene_ontology_full.csv
migraine_gene_ontology_genes.csv
migraine_ALL_SOURCES_GENES.csv
migraine_SOURCES_SUMMARY.csv
```

Source-specific full tables intentionally do **not** use one universal schema because the underlying databases expose different evidence types and metadata. Genes-only files provide a simplified representation for downstream harmonisation.

Examples of source-specific fields include:

- PubMed: `Gene`, `Entrez_ID`, `PMID_Count`, `Supporting_PMIDs`, `Phenotype`, `Source`
- ClinVar: `Gene`, `ClinVar_VariationID`, `ClinVar_Accession`, `Matched_Condition`, `Clinical_Significance`, `Source`
- KEGG: `Phenotype`, `Pathway_ID`, `Pathway_Name`, `KEGG_Gene_ID`, `Gene`, `Source`
- GWAS Catalog: gene symbol, variant identifier, genomic position and mapping context

---

## Gene-symbol harmonisation

The downstream Python workflow validates retrieved identifiers against the NCBI human gene reference.

- current official human gene symbols are retained directly;
- recognised aliases and historical symbols are resolved to current official symbols;
- unresolved records are retained separately rather than silently discarded; and
- contributing-source information is preserved for downstream comparison and ranking.

In the final manuscript benchmark:

- **71,556** combined input symbols were submitted to harmonisation;
- **66,549** were retained after validation (**93.0%**);
- **63,860** matched current official symbols directly; and
- **2,689** were rescued through synonym resolution.

---

## Evidence classes

PhenoGeneR distinguishes evidence types because retrieved genes do not all represent the same biological relationship.

The manuscript analysis uses the following evidence classes:

- curated clinical/phenotype evidence;
- human genetic association evidence;
- integrated disease association evidence;
- literature-derived evidence;
- expression/tissue context;
- pathway or functional annotation;
- network expansion; and
- protein annotation.

This distinction is used to separate broad retrieval from more restrictive direct-evidence prioritisation.

---

## Evidence-aware ranking

The manuscript evaluates the following transparent aggregation strategies:

| Method | Description |
|---|---|
| `Union` | Unranked union of ranking-eligible genes |
| `Source_Count_All` | Number of supporting ranking-eligible sources |
| `Source_Size_Normalized` | Downweights very large source-specific gene sets |
| `Normalized_RRF` | Normalised reciprocal-rank fusion |
| `Direct_Evidence_Count` | Source-count ranking restricted to direct-evidence sources |
| `Direct_Normalized_RRF` | RRF restricted to direct-evidence sources |

These rankings are **prioritisation baselines**, not calibrated probabilities of causality.

---

## Manuscript benchmark

The locked benchmark contains 13 phenotype, disease or trait terms:

```text
asthma
blood pressure medication
body mass index
cholesterol lowering medication
depression
gastro-oesophageal reflux
allergic rhinitis
high cholesterol
hypertension
hypothyroidism
irritable bowel syndrome
migraine
osteoarthritis
```

### Key benchmark results

Across the 12 evaluated sources:

- **84,045** source-level gene retrievals were generated;
- **99/156** source-query combinations returned one or more genes (**63.5%**);
- the mean non-zero pairwise Jaccard similarity was **0.027**, demonstrating substantial source heterogeneity;
- broad retrieval and candidate prioritisation were analysed separately.

Low overlap is interpreted as **source heterogeneity**, not as proof of biological complementarity.

---

## Leakage-controlled benchmark

To reduce circularity, genes supported by at least two of HPO, ClinVar and OMIM were used to construct an internal curated-source consensus reference.

For evaluation:

- HPO was excluded;
- ClinVar was excluded;
- OMIM was excluded; and
- STRING-DB was additionally excluded because its evaluated implementation performs seed-based interaction expansion from curated genes.

Ten benchmark terms had a usable consensus reference, comprising **265 reference genes**.

The remaining eligible evidence sources recovered:

- **174/265** reference genes;
- micro-averaged recovery: **65.7%**;
- macro-averaged recovery: **72.1%**;
- macro Recall@10: **19.9%**;
- macro MRR: **0.384**; and
- macro nDCG@100: **0.270**.

Detailed per-term results are provided in the manuscript Supplementary Data.

---

## Comparator evaluation

### Phen2Gene

Phen2Gene was evaluated only when an exact HPO mapping was available from the benchmark HPO output.

- applicable benchmark terms: **7/13**
- successful results among applicable terms: **7/7**
- median returned genes: **14,621**

Terms without an exact HPO mapping were treated as **not applicable**, not as failures.

### Phenolyzer

Phenolyzer was evaluated using the corresponding phenotype, disease or trait text term.

- successful result production: **11/13**
- median returned genes: **14,967**

The comparison focuses on coverage, candidate-set size, ranking behaviour and interpretability rather than claiming that one method is universally superior.

---

## Limited external ClinGen benchmark

ClinGen Gene-Disease Validity records were used as an external source not queried by the PhenoGeneR retrieval pipeline.

After manual review, six related or subtype associations were retained across two benchmark terms:

### High cholesterol

- `LDLR`
- `LDLRAP1`
- `PCSK9`
- `APOB`

### Migraine

- `SCN1A`
- `ATP1A2`

This is intentionally reported as a **limited external benchmark**, not as validation across all 13 benchmark terms.

In the final analysis:

- `Direct_Normalized_RRF` recovered all six reference genes overall;
- Phenolyzer recovered all six reference genes overall; and
- Phen2Gene was applicable to migraine and recovered both migraine reference genes overall.

Full rank-based metrics are reported in the manuscript Supplementary Data.

---

## Reproducing the manuscript analyses

The final Bioinformatics submission is based on the following analysis stages:

| Analysis | Purpose |
|---|---|
| Analysis 1 | Source/query retrieval coverage |
| Analysis 2 | NCBI symbol validation and synonym rescue |
| Analysis 3 | Cross-source heterogeneity and overlap |
| Analysis 8 | Leakage-controlled recovery benchmark |
| Analysis 10 | Evidence-aware ranking strategies |
| Analysis 11 | Phen2Gene comparator |
| Analysis 12 | Multi-method comparator summary |
| Analysis 13 | Limited external ClinGen benchmark |

Machine-readable analysis outputs are stored under:

```text
AllAnalysisGene/
```

Legacy exploratory analyses may remain in the repository for transparency, but they are **not used as the principal validation evidence in the final manuscript**.

---

## Project structure

A simplified repository layout is:

```text
PhenoGeneR/
├── README.md
├── LICENSE
├── requirements.R
├── requirements.txt
├── environment.yml
├── download_genes.R
├── download_genes_analysis.py
│
├── clinvar.R
├── gene_ontology.R
├── gtex.R
├── gwasrapidd.R
├── hpo.R
├── kegg.R
├── omim.R
├── opentargets.R
├── pubmed.R
├── reactome_pathways.R
├── string_db.R
├── uniprot.R
├── disgenet.R                 # optional/development; not in final benchmark
│
├── Analysis_1_Coverage.py
├── Analysis_2_Validation.py
├── Analysis_3_Overlap.py
├── Analysis_8_*.py
├── Analysis_10_*.py
├── Analysis_11_*.py
├── Analysis_12_*.py
├── Analysis_13_*.py
│
├── AllPackagesGenes/
└── AllAnalysisGene/
```

---

## Reproducibility

For manuscript reproduction, keep the following together in the frozen release:

- exact source code used for the submission;
- R and Python dependency specifications;
- the 13 benchmark terms;
- machine-readable benchmark outputs;
- ranking/comparator outputs;
- test/example data;
- release manifest; and
- checksums for archived benchmark files.

The frozen Bioinformatics submission release should be archived in a persistent repository.

### Persistent archive

**Zenodo DOI:** `TO BE ADDED BEFORE SUBMISSION`

After the DOI is minted, replace the placeholder above and use the same DOI in the manuscript Availability and Data Availability statements.

---

## Important interpretation notes

PhenoGeneR should be interpreted as a **source-aware retrieval, harmonisation and evidence-comparison framework**.

The workflow does **not** claim that:

- every retrieved gene is directly associated with the input phenotype;
- pathway, expression, literature and network evidence are equivalent to curated clinical associations;
- low source overlap proves biological complementarity; or
- its ranking scores are calibrated probabilities of disease causality.

Users should inspect source provenance, evidence class and the biological context of retrieved genes before downstream interpretation.

---

## Troubleshooting

### Missing R packages

```bash
Rscript requirements.R
```

For Bioconductor packages, install `BiocManager` if needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
```

### Empty source results

An empty result can reflect:

- no matching records for the query term;
- source-specific terminology;
- a temporary API/service issue; or
- source-specific query limitations.

Check the source-level output and execution log before interpreting an empty result biologically.

### Network/API problems

If a source is temporarily unavailable, rerun the affected module rather than assuming that no biological evidence exists.

### Optional DisGeNET access

```bash
export DISGENET_API_KEY="YOUR_KEY"
Rscript disgenet.R migraine
```

Never store API credentials directly in tracked source files.

---

## License

PhenoGeneR is released under the [MIT License](LICENSE).

---

## Citation

If you use PhenoGeneR, please cite the accompanying manuscript:

> Muneeb M, Ascher DB. **PhenoGeneR: source-aware retrieval and harmonisation of phenotype-associated gene evidence.** *Bioinformatics*. Application Note. [Citation to be updated after publication.]

---

## Primary resource citations

Users should also cite the underlying resources relevant to their analysis:

- ClinVar
- Human Phenotype Ontology
- OMIM
- GWAS Catalog
- Open Targets
- PubMed / PubTator3
- GTEx
- KEGG
- Reactome
- STRING
- UniProt
- Gene Ontology

The manuscript and Supplementary Data provide the corresponding literature references.

---

## Authors

**Muhammad Muneeb**  
School of Chemistry and Molecular Biosciences  
The University of Queensland  
Brisbane, Australia  
Email: [m.muneeb@uq.edu.au](mailto:m.muneeb@uq.edu.au)

**David B. Ascher**  
The University of Queensland  
Baker Heart and Diabetes Institute  
Email: [d.ascher@uq.edu.au](mailto:d.ascher@uq.edu.au)

---

## Support

For reproducible bug reports and feature requests, use the repository issue tracker:

https://github.com/MuhammadMuneeb007/PhenoGeneR/issues

When reporting an issue, please include:

- operating system;
- R and Python versions;
- command executed;
- complete error message; and
- a minimal reproducible example where possible.

---

## Submission release

**Software name:** PhenoGeneR  
**Licence:** MIT  
**Repository:** https://github.com/MuhammadMuneeb007/PhenoGeneR  

