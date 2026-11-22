"""Auditable DGIdb lookup for theca stress-pathway genes and active TFs.

The output is interaction evidence, not an inferred PCOS treatment list.
"""
import csv
import json
import urllib.request
from pathlib import Path

out = Path('validation/biological_case_Theca_7942968')
de = {r['gene']: r for r in csv.DictReader((out / 'fine5_PCOS_vs_control_DE.csv').open())}
tf = {r['TF']: r for r in csv.DictReader((out / 'fine5_TRRUST_TF_PCOS_vs_control.csv').open())}
ora = {r['ID']: r for r in csv.DictReader((out / 'fine5_Reactome_DEG_ORA_up.csv').open())}
stress_pathways = [
    'REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS',
    'REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES',
]
stress_genes = set()
for pathway in stress_pathways:
    if float(ora[pathway]['p.adjust']) >= .05:
        raise RuntimeError(f'Pathway {pathway} did not pass FDR')
    stress_genes.update(ora[pathway]['geneID'].split('/'))
active_tf = {name for name, row in tf.items()
             if float(row['adj.P.Val']) < .05 and name in ('NFE2L2', 'RELA')}
targets = sorted(stress_genes | active_tf)
endpoint = 'https://dgidb.org/api/graphql'
nodes = []
for i in range(0, len(targets), 5):
    batch = targets[i:i+5]
    names = ', '.join(json.dumps(x) for x in batch)
    query = ('{ genes(names: [' + names + ']) { nodes { name interactions {'
             ' drug { name conceptId approved } interactionScore'
             ' interactionTypes { type directionality }'
             ' sources { sourceDbName } publications { pmid }'
             ' } } } }')
    req = urllib.request.Request(endpoint,
        data=json.dumps({'query': query}).encode(),method='POST',
        headers={'Content-Type':'application/json','Accept':'application/json',
                 'User-Agent':'PCOS-CAM academic analysis'})
    with urllib.request.urlopen(req,timeout=60) as response:
        body = json.load(response)
    if body.get('errors'):
        raise RuntimeError(body['errors'])
    nodes.extend(body['data']['genes']['nodes'])
(out / 'DGIdb_stress_target_raw.json').write_text(
    json.dumps({'source': endpoint, 'targets': targets, 'nodes': nodes},indent=2),
    encoding='utf-8')
rows = []
for node in nodes:
    gene = node['name']
    for record in node.get('interactions') or []:
        drug = record.get('drug') or {}
        rows.append({
            'gene': gene,
            'target_evidence': ';'.join(filter(None,[
                'DE_FDR' if gene in stress_genes and gene in de and
                            float(de[gene]['FDR'])<.05 else '',
                'TF_activity_FDR' if gene in active_tf else '',
                'Reactome_stress_pathway' if gene in stress_genes else ''])),
            'PCOS_vs_control_logFC': de.get(gene,{}).get('logFC',''),
            'gene_DE_FDR': de.get(gene,{}).get('FDR',''),
            'TF_activity_FDR': tf.get(gene,{}).get('adj.P.Val',''),
            'drug': drug.get('name',''),
            'drug_concept_id': drug.get('conceptId',''),
            'DGIdb_approved': drug.get('approved',''),
            'interaction_score': record.get('interactionScore',''),
            'interaction_types': ';'.join(sorted({x.get('type') or '' for x in
                                                  record.get('interactionTypes') or []})),
            'directionality': ';'.join(sorted({x.get('directionality') or '' for x in
                                               record.get('interactionTypes') or []})),
            'source_databases': ';'.join(sorted({x.get('sourceDbName') or '' for x in
                                                 record.get('sources') or []})),
            'pubmed_ids': ';'.join(sorted({str(x.get('pmid')) for x in
                                          record.get('publications') or [] if x.get('pmid')})),
        })
rows.sort(key=lambda r: (r['gene'],-(float(r['interaction_score'])
                                    if r['interaction_score'] else 0)))
columns = list(rows[0]) if rows else ['gene','target_evidence','drug']
with (out / 'DGIdb_stress_target_interactions.csv').open('w',newline='') as handle:
    writer = csv.DictWriter(handle,fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
print('Targets queried:',len(targets),'interactions:',len(rows),
      'approved-drug records:',sum(r['DGIdb_approved'] is True for r in rows))
