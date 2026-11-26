#!/usr/bin/env bash
set -euo pipefail
mkdir -p work out
for t in patient sample; do
  curl -fLsS --retry 2 "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric/data_clinical_$t.txt" -o "work/met_$t.tsv"
done
for g in 898 9134; do
  curl -fLsS --retry 2 "https://www.cbioportal.org/api/molecular-profiles/brca_metabric_mrna/molecular-data?entrezGeneId=$g&sampleListId=brca_metabric_all" -o "work/met$g.json"
done
scan=GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz
curl -fLsS --retry 2 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/$scan" -o "work/$scan"
curl -fLsS --retry 2 'https://raw.githubusercontent.com/12379Monty/GSE96058/master/data/GSE96058_sampDesc.rda' -o work/scan.rda
gzip -dc "work/$scan" | awk -F, 'NR==1 || $1=="\"CCNE1\"" || $1=="\"CCNE2\""{print}' >work/scanexpr.csv
test "$(wc -l < work/scanexpr.csv)" -eq 3
unzip -oq upload/CCNE.zip -d work
chmod -R u+rwX work/CCNE
for f in work/CCNE/*/*.zip; do
  unzip -oq "$f" -d "$(dirname "$f")"
done
