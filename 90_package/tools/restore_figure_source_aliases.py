"""Optionally restore exact duplicate source_data paths expected by old manifests."""
from pathlib import Path
import argparse,csv,hashlib,shutil
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--evidence-root',type=Path,default=Path(__file__).resolve().parents[2]/'05_existing_results')
    ap.add_argument('--restore',action='store_true',help='Materialize copies; without this flag, verify only')
    args=ap.parse_args();root=args.evidence_root.resolve();count=0
    with (root/'source_data_index.tsv').open(encoding='utf-8-sig',newline='') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            src=(root/row['canonical_table']).resolve();dst=(root/row['historical_figure_source_path']).resolve()
            if not src.is_relative_to(root) or not dst.is_relative_to(root):raise ValueError('Path outside evidence root')
            if hashlib.sha256(src.read_bytes()).hexdigest()!=row['sha256']:raise ValueError(src)
            if args.restore:
                if dst.exists() and dst.read_bytes()!=src.read_bytes():raise ValueError('Refusing overwrite '+str(dst))
                dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
            count+=1
    print(f'Verified {count} canonical sources; restored={args.restore}')
if __name__=='__main__':main()
