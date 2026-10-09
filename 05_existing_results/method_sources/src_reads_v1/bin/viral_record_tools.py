#!/usr/bin/env python3
"""Record locked Pyrodigal-gv container and packaged model file inventory."""
import hashlib
import json
from pathlib import Path
import subprocess


def main():
    root = Path(__file__).resolve().parents[1]
    image = json.loads((root / 'bootstrap/reads_receipts/viral-gene.image.json').read_text())[0]
    code = """
import importlib.metadata,json,hashlib
from pathlib import Path
import pyrodigal_gv
d=importlib.metadata.distribution('pyrodigal-gv')
files=[]
for f in d.files:
 p=Path(d.locate_file(f))
 if p.is_file() and ('pyrodigal_gv' in str(f)) and not str(f).endswith('.pyc'):
  files.append(dict(path=str(f),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
print(json.dumps(dict(pyrodigal_gv=d.version,pyrodigal=importlib.metadata.version('pyrodigal'),files=sorted(files,key=lambda x:x['path']))))
"""
    result = subprocess.check_output(['docker', 'run', '--rm', image['Id'], 'python', '-c', code], text=True)
    inventory = json.loads(result)
    doc = dict(schema_version=1, container=image['Id'], local_tag='crc-phage/viral-gene:locked',
               tools={k: inventory[k] for k in ['pyrodigal_gv', 'pyrodigal']},
               model_package_inventory=inventory['files'],
               pixi_lock_sha256=hashlib.sha256((root / 'envs/viral-gene/pixi.lock').read_bytes()).hexdigest(),
               mode='ViralGeneFinder(meta=True, viral_only=False)',
               sources=['https://github.com/althonos/pyrodigal-gv', 'https://pyrodigal.readthedocs.io/en/stable/api/genes.html'],
               scientific_calibration=False)
    doc['model_inventory_sha256'] = hashlib.sha256(json.dumps(inventory, sort_keys=True).encode()).hexdigest()
    target = root / 'assets/viral_gene_tools.json'
    content = json.dumps(doc, sort_keys=True, indent=2) + '\n'
    if not target.exists() or target.read_text() != content:
        target.write_text(content)
    print(image['Id'])


if __name__ == '__main__':
    main()
