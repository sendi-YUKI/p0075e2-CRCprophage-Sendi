#!/usr/bin/env python3
"""Materialize complete result bundles; rewrite only known filesystem roots, never biological IDs."""
import argparse,json,shutil,sys
from pathlib import Path
from reads_io import dump
def rewrite(obj,roots):
    if isinstance(obj,dict):return {rewrite(k,roots):rewrite(v,roots) for k,v in obj.items()}
    if isinstance(obj,list):return [rewrite(x,roots) for x in obj]
    if isinstance(obj,str):
        for src,dst in sorted(roots.items(),key=lambda x:-len(x[0])):
            if obj==src or obj.startswith(src+'/'):return dst+obj[len(src):]
    return obj
def main():
    if len(sys.argv)==3 and not sys.argv[1].startswith('--'):
        root=Path(sys.argv[1]).resolve();target=str(Path(sys.argv[2]).resolve())
        for p in root.rglob('*.json'):
            try:dump(p,rewrite(json.loads(p.read_text()),{str(root):target}))
            except (UnicodeError,json.JSONDecodeError):pass
        return
    p=argparse.ArgumentParser();p.add_argument('--assembly-dirs',nargs='*',required=True);p.add_argument('--outdir',required=True);p.add_argument('--publish-root',required=True);a=p.parse_args()
    if not Path(a.publish_root).is_absolute():raise ValueError('publish-root must be absolute; resolve relative output against Nextflow launchDir')
    out=Path(a.outdir).resolve();out.mkdir(parents=True,exist_ok=True);roots={};files=[]
    for d in a.assembly_dirs:
        d=Path(d).resolve();u=json.load(open(d/'assembly_unit.json'));target=out/'units'/u['unit_id'];target.parent.mkdir(exist_ok=True)
        shutil.copytree(d,target,dirs_exist_ok=True,ignore=shutil.ignore_patterns("readback.sam"));roots[str(d)]=str(Path(a.publish_root)/'units'/u['unit_id'])
        # An annotation bundle points back to the complete original assembly; retain that too.
        if u.get('assembly_fasta') and not Path(u['assembly_fasta']).is_relative_to(d):
            source=Path(u['assembly_fasta']).parent;shutil.copytree(source,target/'assembly_source',dirs_exist_ok=True);roots[str(source)]=str(Path(a.publish_root)/'units'/u['unit_id']/'assembly_source')
        files.append(target/'assembly_unit.json')
    for file in files:dump(file,rewrite(json.load(open(file)),roots))
    from reads_assembly import collect
    collect(files,out)
if __name__=='__main__':main()
