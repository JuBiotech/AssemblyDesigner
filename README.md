# AssemblyDesigner — design & simulation for combinatorial liquid handling and in-silico plasmid design

*Supports MoClo (Golden Gate), Gibson, and 3‑G cloning protocols. Generates worklists for liquid-handling robots and ships an in‑silico assembly pipeline (“Plasmidio”).*

- ✅ Robots: **Opentrons (OT-2 / Flex)**, **Tecan Fluent/EVO**
- 📦 GitHub: [JuBiotech/AssemblyDesigner](https://github.com/JuBiotech/AssemblyDesigner)

## 🚀 Quick start

```bash
git clone https://github.com/JuBiotech/AssemblyDesigner.git
cd AssemblyDesigner
conda create -n adesigner python=3.11 -y
conda activate adesigner
pip install -e ".[dev]"
```

For the optional in-silico stack, install:

```bash
pip install -e ".[insilico]"
```

---

## 📖 Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Example Notebooks](#example-notebooks)
- [How It Works](#how-it-works)
- [Plasmidio (in‑silico) Quickstart](#plasmidio-in-silico-quickstart)
- [3‑G Pipeline](#3-g-pipeline)
- [Installation](#installation)
- [Contributing](#contributing)


---

## 🧩 Overview
`AssemblyDesigner` was developed for **biofoundries** to rapidly generate **liquid-handling worklists** for high-throughput plasmid assembly.
The goal is to enable efficient construction of plasmids for **protein overexpression** and **genomic modifications** in a scalable and automated way.

The toolkit provides:
- **Simulation of liquid-handling routines** → validate volumes, transfers, and mixing steps before execution
- **Automated worklist generation** for plasmid assembly on **Opentrons** and **Tecan** platforms
- A complete **3-G assembly pipeline** (Golden Gate → PCR → Gibson), including


In addition, the package includes a powerful **in-silico module (“Plasmidio”)** that offers:
- **In-silico plasmid assembly** and clean **GenBank outputs**
- **High-throughput sequencing (HTS) analysis** and automated **result plotting**
- **Feature de-duplication and clean-up** for reliable construct annotation
- Simulation of **Golden Gate, Gibson, and PCR assemblies**
- **Provenance tracking and visualization** of assembly histories
- An **automated S1 documentation generator** for standardized experimental records

🎯 **Target audience & goal:**
This module is designed for **scientists with little programming experience** who primarily work in **Jupyter notebooks**.
It substitutes **time-consuming in-silico tasks** that would otherwise distract from laboratory work, providing an accessible, notebook-friendly interface to accelerate **design, validation, and documentation** in synthetic biology.

---

## ✨ Features
- **Worklist generation** for Opentrons & Tecan platforms
- **Simulation** of liquid-handling steps (volumes, mixes, transfers)
- **Plasmidio** (in‑silico):
  - Build assembly reports (`*_report.zip`) and extract **.gb/.gbk** files
  - Feature de‑duplication and clean‑up
  - **3‑G pipeline** (Golden Gate → PCR → Gibson) via `run_3g_batch_safe`
  - Provenance DAGs with `AssemblyHistory`
- Notebook‑friendly, CI‑ready (pytest, pre‑commit)

---

## 📓 Example Notebooks
This repository includes several Jupyter notebooks demonstrating the workflows:

- **GelSim_example.ipynb** → Example simulation of pipetting workflows
- **S1_docu.ipynb** → Generates S1‑style documentation of assemblies
- **Analyzing_SequencingData.ipynb** → Analyze high‑throughput sequencing data for construct verification
- **3G Assembly.ipynb** → End‑to‑end 3‑G pipeline (Golden Gate → PCR → Gibson)
- **AGGA_AutomatedGoldenGateAssembly_.ipynb** → Automated Golden Gate Assembly worklists
- **MoClo WL generator_V8.ipynb** → MoClo worklist generation pipeline
- **PCR_Designer.ipynb** → Primer and PCR fragment design tool

---

## ⚙️ How It Works
1. **Fill Excel template** with design parameters (sources, volumes, mappings).
2. **Simulate** the pipetting steps to validate your plan.
3. **Generate worklists** for robot execution.
4. *(Optional)* **Run Plasmidio** to assemble plasmids in silico and export cleaned GenBank files; analyze sequencing reads; generate documentation.

> **Note:** This repository generates **worklists** (Excel/CSV). Execution scripts for robots depend on your hardware.

---

## 🧬 Plasmidio (in‑silico) Quickstart

```python
from pathlib import Path
from assembly_designer.plasmidio import (
    generate_assembly_reports,
    organize_assembly_reports,
    remove_near_duplicate_features,
)

BASE_DIR = Path.cwd()
folders = [
    BASE_DIR / "Promoter_parts",
    BASE_DIR / "RBS_parts",
    BASE_DIR / "Gene_of_interest_parts",
    BASE_DIR / "Terminator_parts",
    BASE_DIR / "Backbone_parts",
]

# 1) Generate *_report.zip into ./reports/
reports = generate_assembly_reports(folders=folders, output_dir="reports")

# 2) Extract .gb/.gbk into ./reports/Assembly (and optionally delete ZIPs)
gb_paths = organize_assembly_reports("reports", reports, delete_zip=True)

# 3) Clean near-duplicate features in place
for gb in gb_paths:
    remove_near_duplicate_features(gb, tolerance=3)
```

```python
from assembly_designer.plasmidio import run_3g_batch_safe, build_histories_for_all_constructs

batch = run_3g_batch_safe(
    category_dirs=folders,
    designs=excel_path,
    category_order=["Promoter","RBS","Gene","Terminator","Backbone"],
    reports_dir=Path("reports"),
)

# Build provenance graphs per ConstructID
histories = build_histories_for_all_constructs(
    constructs_df=batch.constructs_df,
    assembly_df=batch.assembly_df,
    tus_df=batch.tus_df,
    vector_rec=batch.vector_rec,
    products=batch.products,
    pcr_results=batch.pcr_results,
    finals=batch.finals,
    plot=True,
    save_png_dir=Path("reports/History"),
)
```


# Contributing & development installation

### 1. Clone the repository
GitHub repository: [JuBiotech/AssemblyDesigner](https://github.com/JuBiotech/AssemblyDesigner)

```bash
git clone https://github.com/JuBiotech/AssemblyDesigner.git
cd AssemblyDesigner
```

This repository is maintained on GitHub at [JuBiotech/AssemblyDesigner](https://github.com/JuBiotech/AssemblyDesigner).

> **⚠️ Windows users:** if cloning or installing fails with a "filename too long" / "path too long" error, enable long path support in Git and re-run the failing step:
> ```bash
> git config core.longpaths true
> ```

### 2. Create and activate a fresh environment (example with conda)
```bash
conda create --name adesigner python=3.11 -y
conda activate adesigner
```

### 3. Install uv
```bash
pip install uv
```

### 4. Perform an editable install of AssemblyDesigner including dev dependencies
```bash
uv pip install -e ".[dev]"
```

### 5. Install and enable pre-commit hooks
```bash
pip install pre-commit
pre-commit install
```

### 6. Install in-silico functions (DnaCauldron + SnapGene)
###    Requires Git (for the snapgene_reader dependency)
```bash
uv pip install "assembly_designer[insilico]"
```

### or:
```bash
uv pip install -e ".[insilico]"
```

# 📦 Supported File Formats

This toolkit sticks to common, well-documented formats so you can slot it into existing workflows. Below is what it **reads** and **writes**, grouped by purpose.

## Sequence & Annotation (input)
- **GenBank**: `.gb`, `.gbk`
  Parsed via Biopython. Circular/linear respected; standard qualifiers (`/label`, `/gene`, `/locus_tag`, `/note`) read.
- **SnapGene**: `.dna`
  Read with `snapgene_reader`. Imports sequence, topology, and features (no SnapGene app required).
- **FASTA**: `.fa`, `.fasta`, `.fna`
  Single or multi-record; no feature annotations.
- **CSV/TSV feature tables** *(optional)*: `.csv`, `.tsv`
  If provided, should include columns like `start,end,strand,type,name` (or equivalent).

## Read Data (input, Plasmidio)
- **FASTQ**: `.fastq`, `.fq` (also `.fastq.gz`, `.fq.gz`)
  Single-end supported. Quality scores are ignored for alignment.

## Design & Metadata (input)
- **Excel**: `.xlsx`
  Used by worklist generators and the 3-G pipeline. Validated sheet/column names; see `examples/`.
- **CSV/TSV**: `.csv`, `.tsv`
  Alternative to Excel for maps, part lists, parameters (UTF-8 expected).

## Worklists & Robot I/O (output)
- **Opentrons / Tecan worklists**: `.csv`, `.gwl`
  Human-readable transfer tables; importable on OT-2/Flex and Fluent/EVO (templates in `examples/`)..

## In-silico Assembly & Sequencing Results (output)
- **GenBank (cleaned/assembled)**: `.gb`, `.gbk`
  Fully annotated plasmids after assembly/cleanup.
- **Plasmid assembly reports**: `*_report.zip`
  Contains HTML/CSV summaries and GenBank exports; unpacked to `reports/Assembly/`.
- **QC tables**: `.csv`
  Per-read alignments, feature coverage, PID metrics.
