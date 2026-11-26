# Comparative survival associations and coexpression pathways of CCNE1 and CCNE2 in breast cancer

Analysis code for the manuscript of the same title, a retrospective study of
TCGA-BRCA, METABRIC and SCAN-B.

## Repository structure

| Path | Analysis | Code author |
|---|---|---|
| `retrieve_external.sh`, `analyze.py` | Main survival analyses, sensitivity analyses, KEGG coexpression summaries, Figures 1–3, Supplementary Tables S1–S10 | Jun-Chae Lee |
| `extra/analyze_clinical_expression.py` | Tumor expression and clinicopathologic associations, Figures 4–6, Supplementary Tables S11–S14 | Jun-Chae Lee |
| `additional/run_all.py` | Proliferation adjustment, paralog contrast, copy number, tumor purity; Tables 4–6, Figure 7, Supplementary Tables S15–S18 | Chaeyun Kim |

## How to reproduce

```bash
pip install -r requirements.txt
```

**Main analyses** – see `README_main_analysis.txt`. The TCGA-BRCA discovery
inputs (OncoLnc-format expression–survival exports, the clinical table and the
LinkedOmics report archive) are not redistributed here; place them in
`upload/` and run `bash retrieve_external.sh`, then `python analyze.py`.
`external_input_sha256.txt` lists checksums of the retrieved public inputs.

**Clinical expression analyses** – see `extra/README.txt`.

**Additional analyses** – see `additional/README.txt`. Fully automatic:
`python additional/run_all.py [workdir]` downloads the public data and writes
results to `workdir/results`.

## Data sources

TCGA-BRCA (OncoLnc; cBioPortal PanCanAtlas for the additional analyses),
METABRIC (cBioPortal) and SCAN-B (NCBI GEO GSE96058). No patient-level data
are included in this repository.

## Citation

If you use this code, please cite the manuscript and the archived release
(see `CITATION.cff`).

## License

MIT – see `LICENSE`.
