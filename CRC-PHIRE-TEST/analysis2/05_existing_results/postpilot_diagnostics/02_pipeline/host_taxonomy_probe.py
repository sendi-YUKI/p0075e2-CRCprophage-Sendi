import pathlib,bz2,pickle,json,csv,hashlib
R=pathlib.Path('/srv/CRC-PHIRE/analysis1');A=R/'14_postpilot_adjustment'
db=pathlib.Path('/srv/CRC-PHIRE/databases/metaphlan/mpa_vJan25_CHOCOPhlAnSGB_202503/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl')
with open(db,'rb') as f:compressed=f.read(3)==b'BZh'
with (bz2.open(db,'rb') if compressed else open(db,'rb')) as f:d=pickle.load(f)
tax=d['taxonomy'];terms=['Bacteroides_fragilis','Bacteroides_hominis','Fusobacterium_animalis','Fusobacterium_nucleatum','Parvimonas_micra','Peptostreptococcus_anaerobius','Faecalibacterium_prausnitzii','Roseburia_intestinalis']
rows=[]
for term in terms:
    matched=[k for k in tax if any(x=='s__'+term for x in k.split('|'))]
    rows.append({'requested_taxon':term,'exact_species_in_database':bool(matched),'SGB_taxonomy_records':len(matched),'matching_clades':';'.join(matched)})
with open(A/'03_host_profile/database_taxonomy_coverage.tsv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
first=next(iter(d['markers'].items()))
(A/'03_host_profile/database_taxonomy_probe.json').write_text(json.dumps({'database':str(db),'database_sha256':hashlib.sha256(db.read_bytes()).hexdigest(),'taxonomy_records':len(tax),'marker_records':len(d['markers']),'marker_example':first,'keys':list(d),'requested_matches':[{k:v for k,v in r.items() if k!='matching_clades'} for r in rows]},indent=2))
print(json.dumps([{k:v for k,v in r.items() if k!='matching_clades'} for r in rows]))
