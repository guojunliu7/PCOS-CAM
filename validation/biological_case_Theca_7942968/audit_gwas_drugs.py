"""Reproducible, conservative target audit from archived public database results.

No database interaction is interpreted as evidence of clinical efficacy.
Run from the PCOS-CAM project root.
"""
import csv
import json
from pathlib import Path

out = Path('validation/biological_case_Theca_7942968')
gwas_file = Path('validation/biological_case_GSE313660/GWAS_Catalog_PCOS_associations_raw.json')


def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, columns):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


de = {r['gene']: r for r in read_csv(out / 'fine5_PCOS_vs_control_DE.csv')}
assoc = json.loads(gwas_file.read_text())['records']
overlap = []
for a in assoc:
    try:
        if float(a['p_value']) > 5e-8:
            continue
    except (KeyError, TypeError, ValueError):
        continue
    for gene in a.get('mapped_genes') or []:
        if gene in de and float(de[gene]['FDR']) < .05:
            overlap.append({
                'gene': gene,
                'GWAS_association_id': a.get('association_id'),
                'GWAS_study': a.get('accession_id'),
                'GWAS_p_value': a.get('p_value'),
                'fine5_logFC': de[gene]['logFC'],
                'fine5_FDR': de[gene]['FDR'],
                'interpretation': 'mapped locus gene, not a proven causal gene',
            })
write_csv(out / 'fine5_GWAS_DE_overlap.csv', overlap,
          ['gene','GWAS_association_id','GWAS_study','GWAS_p_value',
           'fine5_logFC','fine5_FDR','interpretation'])

drugs = read_csv(out / 'DGIdb_stress_target_interactions.csv')
approved = [r for r in drugs if r['DGIdb_approved'] == 'True']
directional = [r for r in approved if r['directionality']]
for r in approved:
    r['PCOS_GWAS_supported_target'] = 'yes' if any(
        x['gene'] == r['gene'] for x in overlap) else 'no'
    r['status'] = 'database interaction only; not a PCOS therapy candidate'
write_csv(out / 'DGIdb_approved_hits_audit.csv', approved,
          list(approved[0].keys()) if approved else ['gene','drug','status'])
write_csv(out / 'DGIdb_directional_approved_hits_audit.csv', directional,
          list(approved[0].keys()) if approved else ['gene','drug','status'])
print('PCOS GWAS / cluster-5 FDR-DE gene overlaps:', len({x['gene'] for x in overlap}))
print('DGIdb total:', len(drugs), 'approved:', len(approved),
      'approved with recorded direction:', len(directional))
print('Approved interactions for a GWAS/DE overlap gene:',
      sum(r['PCOS_GWAS_supported_target'] == 'yes' for r in approved))
