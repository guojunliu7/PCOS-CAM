"""Screen GWAS Catalog PCOS associations against exploratory macrophage results.

Mapped genes are variant-overlap/proximity annotations, not causal target claims.
The script intentionally does not infer a drug candidate from nominal RNA signals.
"""
import csv
import json
import urllib.parse
import urllib.request
from pathlib import Path

out = Path('validation/biological_case_GSE313660')
endpoint = 'https://www.ebi.ac.uk/gwas/rest/api/v2/associations'
query = urllib.parse.urlencode({'efo_trait': 'polycystic ovary syndrome', 'size': 500})
url = endpoint + '?' + query
req = urllib.request.Request(url, headers={'Accept': 'application/json',
                                            'User-Agent': 'PCOS-CAM academic analysis'})
with urllib.request.urlopen(req, timeout=60) as response:
    data = json.load(response)
records = data.get('_embedded', {}).get('associations', [])
total = data['page']['totalElements']
if len(records) != total:
    raise RuntimeError(f'Expected {total} associations; retrieved {len(records)}')
(out / 'GWAS_Catalog_PCOS_associations_raw.json').write_text(
    json.dumps({'source': url, 'records': records}, indent=2), encoding='utf-8')
de = {}
with (out / 'macrophage_PCOS_vs_control_DE.csv').open(newline='') as handle:
    for row in csv.DictReader(handle):
        de[row['gene']] = row
rows = []
for association in records:
    p = association.get('p_value')
    if p is None or float(p) > 5e-8:
        continue
    for gene in association.get('mapped_genes') or []:
        if gene not in de:
            continue
        row = de[gene]
        rows.append({
            'gene': gene,
            'GWAS_association_id': association['association_id'],
            'GWAS_study': association.get('accession_id'),
            'GWAS_pubmed_id': association.get('pubmed_id'),
            'GWAS_p_value': p,
            'GWAS_location': ';'.join(association.get('locations') or []),
            'GWAS_gene_mapping': 'mapped_genes; locus annotation, not causal proof',
            'macrophage_logFC_PCOS_vs_control': row['logFC'],
            'macrophage_PValue': row['PValue'],
            'macrophage_FDR': row['FDR'],
            'passes_macrophage_FDR_0_05': float(row['FDR']) < 0.05,
        })
columns = ['gene', 'GWAS_association_id', 'GWAS_study', 'GWAS_pubmed_id',
           'GWAS_p_value', 'GWAS_location', 'GWAS_gene_mapping',
           'macrophage_logFC_PCOS_vs_control', 'macrophage_PValue',
           'macrophage_FDR', 'passes_macrophage_FDR_0_05']
with (out / 'GWAS_macrophage_gene_overlap.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
print('GWAS associations:', len(records), 'overlap rows:', len(rows),
      'distinct genes:', len({r['gene'] for r in rows}),
      'DE FDR-supported:', sum(r['passes_macrophage_FDR_0_05'] for r in rows))
