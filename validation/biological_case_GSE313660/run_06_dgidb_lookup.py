"""Record DGIdb interactions for PCOS GWAS-mapped genes expressed in macrophages.

This is an audit of database interactions, not a nominated PCOS drug list.
Only genes passing patient-level DE FDR would be eligible for nomination.
"""
import csv
import json
import urllib.request
from pathlib import Path

out = Path('validation/biological_case_GSE313660')
overlap = list(csv.DictReader((out / 'GWAS_macrophage_gene_overlap.csv').open()))
genes = sorted({r['gene'] for r in overlap})
de = {r['gene']: r for r in overlap}
endpoint = 'https://dgidb.org/api/graphql'
all_nodes = []
for i in range(0, len(genes), 5):
    batch = genes[i:i + 5]
    names = ', '.join(json.dumps(g) for g in batch)
    query = ('{ genes(names: [' + names + ']) { nodes { name interactions {'
             ' drug { name conceptId } interactionScore'
             ' interactionTypes { type directionality }'
             ' sources { sourceDbName } publications { pmid }'
             ' } } } }')
    req = urllib.request.Request(
        endpoint, data=json.dumps({'query': query}).encode(), method='POST',
        headers={'Content-Type': 'application/json', 'Accept': 'application/json',
                 'User-Agent': 'PCOS-CAM academic analysis'})
    with urllib.request.urlopen(req, timeout=60) as response:
        body = json.load(response)
    if body.get('errors'):
        raise RuntimeError(body['errors'])
    all_nodes.extend(body['data']['genes']['nodes'])
(out / 'DGIdb_GWAS_genes_raw.json').write_text(
    json.dumps({'source': endpoint, 'genes': genes, 'nodes': all_nodes}, indent=2),
    encoding='utf-8')
rows = []
for node in all_nodes:
    gene = node['name']
    if gene not in de:
        continue
    for item in node.get('interactions') or []:
        drug = item.get('drug') or {}
        rows.append({
            'gene': gene, 'drug': drug.get('name'),
            'drug_concept_id': drug.get('conceptId'),
            'interaction_score': item.get('interactionScore'),
            'interaction_types': ';'.join(sorted({x.get('type') or '' for x in
                                                  item.get('interactionTypes') or []})),
            'directionality': ';'.join(sorted({x.get('directionality') or '' for x in
                                               item.get('interactionTypes') or []})),
            'source_databases': ';'.join(sorted({x.get('sourceDbName') or '' for x in
                                                 item.get('sources') or []})),
            'pubmed_ids': ';'.join(sorted({str(x.get('pmid')) for x in
                                          item.get('publications') or [] if x.get('pmid')})),
            'macrophage_logFC': de[gene]['macrophage_logFC_PCOS_vs_control'],
            'macrophage_FDR': de[gene]['macrophage_FDR'],
            'eligible_for_prioritized_drug_list': float(de[gene]['macrophage_FDR']) < .05,
        })
rows.sort(key=lambda r: (r['gene'], -(float(r['interaction_score'])
                                     if r['interaction_score'] else 0)))
columns = list(rows[0]) if rows else [
    'gene', 'drug', 'drug_concept_id', 'interaction_score', 'interaction_types',
    'directionality', 'source_databases', 'pubmed_ids', 'macrophage_logFC',
    'macrophage_FDR', 'eligible_for_prioritized_drug_list']
with (out / 'DGIdb_interactions_exploratory.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
with (out / 'DGIdb_prioritized_drugs.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(r for r in rows if r['eligible_for_prioritized_drug_list'])
print('GWAS-mapped genes queried:', len(genes),
      'DGIdb interactions:', len(rows),
      'prioritized drug interactions:',
      sum(r['eligible_for_prioritized_drug_list'] for r in rows))
