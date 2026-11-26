CCNE1 and CCNE2 in breast cancer
Reproducibility materials

Figures 1–3 are provided at 330 dpi. All tables contain aggregate results.
cox.csv includes model-level N, events, HR, 95% CI, P, expression-term
proportional hazards tests, and apparent concordance where applicable. The
external_q column contains BH-adjusted P values across the four main external
tests. Other q values are within cohort and model, and the sensitivity
results should be interpreted as exploratory. model_shape.csv contains
spline-versus-linear likelihood ratio tests and changes in apparent C.
kegg_all_pathways.csv combines the top-pathway and remaining background
records from the supplied HTML reports. kegg_exports.csv preserves the
displayed records and leading-edge annotations. kegg_overlap.csv compares
the pathways with displayed leading-edge lists for both genes.

To reproduce, put the four supplied analysis inputs (three CSVs and CCNE.zip)
in upload/ in a clean working directory. Run bash retrieve_external.sh to
retrieve public cohorts and expand the pathway exports. Then run analyze.py
from the same working directory. The original run used Python 3.12 with
lifelines 0.30.3, rdata 1.1.0, pandas 2.3.3, scipy 1.18.1, numpy 2.5.4,
and matplotlib 3.11.2. Public cohort data retrieved on 2026-11-05. Future
source revisions can affect counts and estimates. The external data
retrieval script requires unzip and about 600 MB of disk space.
