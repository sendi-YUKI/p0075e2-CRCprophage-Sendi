"""Redraw a real pilot result from delivered TSV only. No server access needed.
Usage: python redraw_replication.py --out-dir /chosen/output
Requires matplotlib. This is a new handoff visualization, not the original figure.
"""
import argparse,csv
from pathlib import Path
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out-dir',type=Path,required=True);a=ap.parse_args()
    root=Path(__file__).resolve().parents[2]
    path=root/'05_existing_results/postpilot_diagnostics/05_replication_activity/replication_candidates.tsv'
    rows=list(csv.DictReader(path.open(encoding='utf-8-sig'),delimiter='\t'))
    # Discover actual ratio field explicitly rather than assuming a guessed table column.
    ratio=next((k for k in ['tool_prophage-host_ratio'] if k in rows[0]),None)
    if ratio is None:raise ValueError('Please set ratio column from table header: '+','.join(rows[0]))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    a.out_dir.mkdir(parents=True,exist_ok=True)
    fig,ax=plt.subplots(figsize=(7,3),layout='constrained')
    ax.scatter([float(r[ratio]) for r in rows],range(len(rows)),s=75,color='#347EA5')
    ax.axvline(2,ls='--',color='#888',lw=1)
    labels=[r.get('donor_id',r.get('sample_id',r.get('biological_sample_id',str(i+1)))) for i,r in enumerate(rows)]
    ax.set_yticks(range(len(rows)),labels);ax.set_xlabel('Prophage / host mean coverage');ax.set_title('Two pilot replication-compatible candidates')
    ax.spines[['top','right']].set_visible(False)
    fig.savefig(a.out_dir/'pilot_replication_candidates.png',dpi=220);fig.savefig(a.out_dir/'pilot_replication_candidates.svg');plt.close(fig)
if __name__=='__main__':main()
