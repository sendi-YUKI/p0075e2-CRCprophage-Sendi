"""Reproducible publication-style pilot figures from verified result tables.

Separate presentation layer: does not change calling, clustering or statistical
thresholds. Production output requires accepted upstream tables. Synthetic layout
previews are explicitly watermarked and kept under a separate root.
"""
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
import argparse,csv,hashlib,html,json,math,os,sys
HERE=Path(__file__).resolve().parent
if (HERE/'figure_preview_dependencies').is_dir():sys.path.insert(0,str(HERE/'figure_preview_dependencies'))
COLORS={'CRC':'#D55E00','control':'#0072B2','reference':'#009E73','viral':'#8064A2','unknown':'#D9DEE3'}
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def table(p):
    with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def export(p,rows):
    rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ['status']
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def selection(tables):
    """Rank displays using pooled occurrence/detection, never CRC contrast."""
    broad={r['m_sequence_id'] for r in tables['matches'] if r['match']!='partial_or_fragment_match'}
    donors=defaultdict(set)
    for r in tables['links']:
        if r['match']!='partial_or_fragment_match':donors[r['m_sequence_id']].add(r['patient_sample'])
    seq={r['viral_sequence_id']:r for r in tables['seq']}
    selected=sorted(broad,key=lambda s:(-len(donors[s]),-int(seq[s]['length']),s))[:16]
    values=defaultdict(list)
    for r in tables['donor']:values[r['reference_id']].append(r)
    def rank(k):
        rows=values[k];det=sum(r['detection_status']=='detected' for r in rows)
        vals=[float(r['fragment_rpkm_all_clean']) for r in rows if r['measurement_status']=='assessed' and r['fragment_rpkm_all_clean']!='']
        return (-det,-sum(vals)/len(vals) if vals else 0,k)
    targets=sorted(values,key=rank)[:24]
    return selected,targets
def panel(ax,letter,title):
    ax.set_title(title,loc='left',fontsize=8,pad=10)
    ax.text(-.10,1.05,letter,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
def empty(ax,text):
    ax.axis('off');ax.text(.5,.5,text,transform=ax.transAxes,ha='center',va='center',wrap=True,color='#435466')
def figure_overview(c):
    plt,np,t=c['plt'],c['np'],c['tables'];groups=c['groups'];report=c['report']
    fig,(ax,bx)=plt.subplots(1,2,figsize=(7.1,3.6),gridspec_kw={'width_ratios':[1,1.3]},layout='constrained')
    ax.axis('off');panel(ax,'a','Study inputs and shared analysis')
    boxes=[(.02,.69,.43,.19,f"{report['reference_genomes']} reference\nbacterial genomes",COLORS['reference']),(.55,.69,.43,.19,f"{len(groups)} stool donors\n{Counter(groups.values())['CRC']} CRC + {Counter(groups.values())['control']} controls",COLORS['control']),(.14,.35,.72,.20,'Shared viral catalog\nSequence links + functions',COLORS['viral']),(.14,.04,.72,.17,'Per-donor DNA measurement\nTechnical / exploratory pilot','#526273')]
    from matplotlib.patches import FancyBboxPatch
    for x,y,w,h,label,color in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.018',facecolor='white',edgecolor=color,lw=1.15,transform=ax.transAxes))
        ax.text(x+w/2,y+h/2,label,ha='center',va='center',transform=ax.transAxes,fontsize=7)
    for start,end in [((.235,.67),(.38,.57)),((.765,.67),(.62,.57)),((.5,.33),(.5,.23))]:ax.annotate('',xy=end,xytext=start,xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#697582',lw=.9))
    panel(bx,'b','Recovered viral instances per donor')
    counts=Counter();provirus=Counter()
    for r in t['occurrences']:
        sid=r['biological_sample_id']
        if sid in groups:
            counts[sid]+=1
            if r['candidate_type']=='provirus_locus':provirus[sid]+=1
    ids=c['donors'];rows=[]
    for y,sid in enumerate(ids):
        bx.plot([provirus[sid],counts[sid]],[y,y],c='#CDD3D8',lw=1.3,zorder=1)
        bx.scatter(counts[sid],y,marker='o',s=30,facecolor=COLORS[groups[sid]],edgecolor='white',lw=.5,zorder=3)
        bx.scatter(provirus[sid],y,marker='|',s=80,color='#172B3A',zorder=4)
        rows.append(dict(sample_id=sid,group=groups[sid],raw_viral_instances=counts[sid],predicted_provirus_loci=provirus[sid]))
    bx.set_yticks(range(len(ids)),ids);bx.invert_yaxis();bx.set_xlim(left=0);bx.set_xlabel('Instances (count)')
    from matplotlib.lines import Line2D
    bx.legend(handles=[Line2D([],[],marker='o',ls='',color=COLORS[g],label=g) for g in ['CRC','control']]+[Line2D([],[],marker='|',ls='',color='#172B3A',label='Provirus subset')],frameon=False,fontsize=6,loc='upper center',bbox_to_anchor=(.5,-.17),ncol=3)
    c['save'](fig,'Fig01_design_and_discovery');c['export']('Fig01_donor_instances',rows)
    return dict(stem='Fig01_design_and_discovery',title='研究输入与逐人发现概览',caption='圆点是每人组装得到的原始病毒实例数，竖线是其中预测provirus位点数。两者不是两个独立集合；原始病毒实例不能全部称为prophage。数量受测序深度和组装能力影响，不能直接解释为疾病差异。',selection='All fixed donors; zeros retained.')
def figure_connections(c):
    plt,np,t=c['plt'],c['np'],c['tables'];chosen=c['selected_m_sequences'];aliases=c['aliases']
    fig,(ax,bx)=plt.subplots(1,2,figsize=(7.1,3.7),gridspec_kw={'width_ratios':[1,1.2]},layout='constrained')
    panel(ax,'a','Direct G-M sequence evidence')
    matches=t['matches'];scatter=None
    for kind,marker in [('partial_or_fragment_match','^'),('broad_bidirectional_match','o'),('identical_oriented_sequence_sha256','s')]:
        rows=[r for r in matches if r['match']==kind]
        if rows:scatter=ax.scatter([float(r['g_coverage_percent']) for r in rows],[float(r['m_coverage_percent']) for r in rows],c=[float(r['ani_percent']) for r in rows],cmap='viridis',vmin=95,vmax=100,marker=marker,s=23,alpha=.7,edgecolors='none',label={'partial_or_fragment_match':'Partial','broad_bidirectional_match':'Broad','identical_oriented_sequence_sha256':'Identical SHA'}[kind])
    ax.axvline(85,color='#87939C',ls='--',lw=.7);ax.axhline(85,color='#87939C',ls='--',lw=.7)
    ax.set(xlim=(0,102),ylim=(0,102),xlabel='G reference sequence covered (%)',ylabel='M stool sequence covered (%)')
    if scatter:
        fig.colorbar(scatter,ax=ax,shrink=.65,pad=.02,label='ANI (%)');ax.legend(frameon=False,loc='lower left',fontsize=6)
    else:ax.text(.45,.5,'No pairs passed the\npredefined initial screen',transform=ax.transAxes,ha='center')
    panel(bx,'b','Reference links for stool sequences')
    edges=defaultdict(set);taxa=set()
    for r in t['links']:
        if r['match']=='partial_or_fragment_match' or r['m_sequence_id'] not in chosen:continue
        taxon=r.get('reference_taxon') or 'Unassigned';taxa.add(taxon);edges[(taxon,r['m_sequence_id'])].add(r['reference_genome'])
    rows=[]
    if chosen and taxa:
        taxa=sorted(taxa);arr=np.zeros((len(taxa),len(chosen)),dtype=int)
        for i,taxon in enumerate(taxa):
            for j,sid in enumerate(chosen):
                arr[i,j]=len(edges[taxon,sid]);rows.append(dict(reference_taxon=taxon,m_sequence_id=sid,alias=aliases[sid],linked_distinct_reference_genomes=int(arr[i,j])))
        im=bx.imshow(arr,cmap='Blues',aspect='auto',vmin=0);bx.set_yticks(range(len(taxa)),taxa,fontsize=6);bx.set_xticks(range(len(chosen)),[aliases[s] for s in chosen],rotation=90,fontsize=6)
        bx.set_xlabel('Stool sequence alias');bar=fig.colorbar(im,ax=bx,shrink=.65,pad=.02,label='Distinct reference genomes (count)')
        from matplotlib.ticker import MaxNLocator
        bar.locator=MaxNLocator(integer=True,nbins=4);bar.update_ticks()
        for i,j in zip(*np.where(arr>0)):bx.text(j,i,str(arr[i,j]),ha='center',va='center',fontsize=5.5,color='white' if arr[i,j]>arr.max()/2 else '#172B3A')
    else:empty(bx,'No broad reference-to-stool\nsequence links available')
    c['save'](fig,'Fig02_sequence_and_reference_links');c['export']('Fig02_direct_pairs',matches);c['export']('Fig02_reference_link_matrix',rows)
    return dict(stem='Fig02_sequence_and_reference_links',title='两端的序列交集与参考菌株连接',caption='左图包含通过既定初筛的直接G–M序列对，可能有重叠点；不是所有未匹配对的分布。右图按参考分类单元和患者病毒序列展示连接，每格按不同参考genome去重计数，不将展开的多对多对应行当独立证据。这里的菌名表示参考端来源，不能证明患者体内宿主。',selection='Up to16 broad-matched M sequences, sorted by pooled distinct-donor count, then length, then ID. No ranking by CRC effect.')
def figure_measurement(c):
    plt,np,t=c['plt'],c['np'],c['tables'];targets=c['targets'];ids=c['donors'];groups=c['groups']
    fig,(ax,bx)=plt.subplots(1,2,figsize=(7.1,5.0),gridspec_kw={'width_ratios':[1.35,1]},layout='constrained')
    panel(ax,'a','Individual donor measurements');panel(bx,'b','Detected repertoire per donor')
    lookup={(r['reference_id'],r['biological_sample_id']):r for r in t['donor']}
    rows=[]
    if targets:
        arr=np.full((len(targets),len(ids)),np.nan)
        for i,target in enumerate(targets):
            for j,sid in enumerate(ids):
                r=lookup.get((target,sid))
                if r and r['measurement_status']=='assessed' and r['fragment_rpkm_all_clean']!='':arr[i,j]=np.log10(1+float(r['fragment_rpkm_all_clean']))
                if r:rows.append(dict(r,display_alias='V'+str(i+1).zfill(2),display_log10_one_plus_rpkm=float(arr[i,j]) if np.isfinite(arr[i,j]) else 'NA'))
        cmap=plt.get_cmap('viridis').copy();cmap.set_bad(COLORS['unknown']);im=ax.imshow(np.ma.masked_invalid(arr),aspect='auto',cmap=cmap,vmin=0)
        ax.set_yticks(range(len(targets)),['V'+str(i+1).zfill(2) for i in range(len(targets))]);ax.set_xticks(range(len(ids)),ids,rotation=90,fontsize=6)
        ax.set_ylabel('vOTU display alias');fig.colorbar(im,ax=ax,orientation='horizontal',pad=.03,shrink=.85,label='log10(1 + fragment RPKM)')
        for j,sid in enumerate(ids):ax.scatter(j,-.88,marker='s',s=22,color=COLORS[groups[sid]],clip_on=False)
    else:empty(ax,'No measurable catalog targets')
    counts=Counter(r['biological_sample_id'] for r in t['donor'] if r['detection_status']=='detected');assessed=Counter(r['biological_sample_id'] for r in t['donor'] if r['detection_status'] in {'detected','not_detected'})
    points=[]
    for x,group in enumerate(['CRC','control']):
        members=[s for s in ids if groups[s]==group];offsets=np.linspace(-.16,.16,len(members)) if len(members)>1 else [0]
        valid=[]
        for offset,sid in zip(offsets,members):
            if assessed[sid]:
                bx.scatter(x+offset,counts[sid],s=35,c=COLORS[group],edgecolor='white',lw=.6,zorder=3);valid.append(counts[sid])
            # Full IDs remain in source data; avoid clutter from point labels.
            points.append(dict(sample_id=sid,group=group,detected_votus=counts[sid] if assessed[sid] else 'NA',detection_evaluable_votus=assessed[sid]))
        if valid:
            median=float(np.median(valid));bx.plot([x-.24,x+.24],[median,median],color='#172B3A',lw=1.4)
        missing=len(members)-len(valid)
        if missing:bx.text(x,.10,f'{missing} donor(s) not assessable',transform=bx.get_xaxis_transform(),ha='center',fontsize=5.5)
    bx.set_xticks([0,1],[f"CRC\n(n={Counter(groups.values())['CRC']})",f"Control\n(n={Counter(groups.values())['control']})"]);bx.set_xlim(-.55,1.55);bx.set_ylim(0,max([1,*counts.values()])*1.32);bx.set_ylabel('Detected vOTUs per donor (count)')
    bx.text(.5,.99,'Each point = one donor\nBlack line = group median',transform=bx.transAxes,ha='center',va='top',fontsize=6.5)
    c['save'](fig,'Fig03_donor_heatmap_and_points');c['export']('Fig03_display_measurements',rows);c['export']('Fig03_all_donor_detection_counts',points)
    return dict(stem='Fig03_donor_heatmap_and_points',title='逐人的病毒谱与检出情况',caption='热图展示按全体供者合并检出人数、合并平均丰度及编号排序的前24个vOTU，未按CRC差异挑选。灰色是未评估，低丰度与零值使用同一连续色阶。右图每个点是一名供者，短横线是中位数；展示全部目录中的检出数量。检出计数受测序和宿主背景影响，图中没有显著性检验或因果含义。',selection='Up to24 pooled-detection-ranked vOTUs; all donors retained. Full unfiltered tables preserved.')
def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=HERE.parent);p.add_argument('--pilot',choices=['pilot6','pilot10'],required=True);p.add_argument('--synthetic-preview',action='store_true');a=p.parse_args()
    root=a.root.resolve();out=root/'12_reports'/a.pilot;catalog=root/'06_viral_catalog'/(a.pilot+'_combined')
    os.environ.setdefault('MPLCONFIGDIR',str(root/'92_tmp/figure_mplconfig'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import font_manager
    font='Arial' if any(f.name=='Arial' for f in font_manager.fontManager.ttflist) else 'DejaVu Sans'
    plt.rcParams.update({'font.family':font,'font.size':7,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.65,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':False,'legend.frameon':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white','text.color':'#172B3A','axes.labelcolor':'#172B3A','savefig.dpi':300})
    acceptance=out/(a.pilot+'_final_acceptance.json');report=read(acceptance)
    assert report['status'] in ['completed_reduced_engineering_pilot','completed_balanced_engineering_pilot']
    assert bool(report.get('synthetic_preview',False))==a.synthetic_preview,'Synthetic and production evidence must never mix'
    if a.synthetic_preview:assert 'figure_style_preview' in str(root)
    else:assert os.environ.get('SLURM_JOB_ID'),'Production figures run through Slurm'
    for rel,h in report['output_sha256'].items():assert sha(out/rel)==h,rel
    frozen=root/'01_manifests'/(a.pilot+'_frozen_assemblies.json');groups=read(frozen)['groups'];assert len(groups)==report['donors']
    files={'matches':out/'G_M_sequence_matches.tsv','links':out/'G_M_occurrence_host_links.tsv','features':out/'intersection_function_evidence.tsv','donor':out/'donor_votu_measurement.tsv','seq':catalog/'catalog/viral_sequence_master.tsv','occurrences':catalog/'catalog/source_occurrences.tsv','genes':catalog/'genes/gene_table.tsv','phrog_hits':catalog/'functions/function_hits.tsv','function_status':catalog/'functions/viral_sequence_function_status.tsv','systems':catalog/'functions/systems.tsv'}
    acceptance_inputs=[]
    if not a.synthetic_preview:
        from pilot10_catalog_acceptance import verified_catalog_completion
        publication=verified_catalog_completion(root,catalog)
        acceptance_inputs=[catalog/'pipeline_info/stage_state.json',catalog/'pipeline_info/terminal_output_manifest.json']
        if publication is not None:acceptance_inputs.append(publication)
        terminal=read(catalog/'pipeline_info/terminal_output_manifest.json')
        for key in ['seq','occurrences']:
            rel=files[key].relative_to(catalog).as_posix();assert rel in terminal,('Unverified figure input',rel)
            assert sha(files[key])==terminal[rel],rel
        for receipt_rel,keys in [('genes/gene_calling_status.json',['genes']),('functions/viral_function_status.json',['phrog_hits','function_status','systems'])]:
            receipt_path=catalog/receipt_rel;assert sha(receipt_path)==terminal[receipt_rel]
            receipt=read(receipt_path);assert receipt['status']=='completed'
            for key in keys:
                name=files[key].name;assert sha(files[key])==receipt['output_files'][name],name
    tables={k:table(v) for k,v in files.items()};figdir=out/'figures';figdir.mkdir(exist_ok=True)
    assert not (figdir/'figure_manifest.json').exists(),'Keep accepted figures immutable; use reviewed new output version to redraw'
    (figdir/'source_data').mkdir(exist_ok=True)
    selected,targets=selection(tables);aliases={s:'M'+str(i+1).zfill(2) for i,s in enumerate(selected)}
    def save(fig,stem):
        if a.synthetic_preview:fig.text(.5,-.10,'SYNTHETIC LAYOUT TEST - NOT STUDY RESULTS',ha='center',fontsize=8,color='#AB3900')
        for ext in ['svg','pdf','png']:fig.savefig(figdir/(stem+'.'+ext),bbox_inches='tight',pad_inches=.10)
        plt.close(fig)
    c=dict(plt=plt,np=np,out=figdir,tables=tables,selected_m_sequences=selected,aliases=aliases,targets=targets,groups=groups,donors=sorted(groups,key=lambda s:(groups[s]!='CRC',s)),report=report,save=save,export=lambda stem,rows:export(figdir/'source_data'/(stem+'.tsv'),rows))
    figures=[figure_overview(c),figure_connections(c),figure_measurement(c)]
    from pilot_function_figures import draw_function_figures
    figures.extend(draw_function_figures(c))
    c['export']('sequence_display_aliases',[dict(alias=v,viral_sequence_id=k) for k,v in aliases.items()])
    c['export']('votu_display_aliases',[dict(alias='V'+str(i+1).zfill(2),reference_id=k) for i,k in enumerate(targets)])
    title=f"CRC-PHIRE · {report['donors']}人先导图文报告";notice='【合成排版预览；不是本项目结果】' if a.synthetic_preview else '真实结果作图 · 小样本技术与描述性先导'
    blocks=[];md=[f'# {title}',notice,'图内英文便于后续论文修改；以下中文说明用于解释。原始分析结果与完整数据表保持不变。']
    for i,f in enumerate(figures,1):
        stem=f['stem'];caption=f['caption'];links=' · '.join(f'<a href="figures/{stem}.{ext}">{ext.upper()}</a>' for ext in ['svg','pdf','png'])
        blocks.append(f'<section><div class="eyebrow">FIGURE {i:02d}</div><h2>{html.escape(f["title"])}</h2><img src="figures/{stem}.svg" alt="{html.escape(f["title"])}"><p>{html.escape(caption)}</p><p class="rule">展示选择：{html.escape(f.get("selection","See figure source data."))}</p><p>{links}</p></section>')
        md.extend([f'## Figure {i}：{f["title"]}',f'![{f["title"]}](figures/{stem}.png)',caption,f'展示规则：{f.get("selection","见源数据")}',f'[可编辑SVG](figures/{stem}.svg) · [PDF](figures/{stem}.pdf)'])
    note='本批样本同时参与发现与测量；图中连接不是已证实的患者宿主，功能同源证据不是实验活性。完整结果与未展示对象保留在原表。'
    page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>body{margin:0;background:#f4f6f7;color:#192d3b;font:16px/1.8 system-ui,"Microsoft YaHei",sans-serif}main{max-width:1120px;margin:auto;padding:38px 24px}header{padding:20px 0 30px;border-bottom:2px solid #192d3b}h1{font-size:30px;line-height:1.3}h2{font-size:23px}.eyebrow{font:12px/1.5 system-ui;letter-spacing:2px;color:#526575}section{background:white;padding:30px;margin:26px 0;border:1px solid #dce3e7;border-radius:8px}img{display:block;width:100%;height:auto;margin:22px auto}.rule{font-size:13px;color:#526575}a{color:#0072b2}footer{font-size:14px}</style><main><header><div class="eyebrow">CRC-PHIRE / PILOT FIGURES</div><h1>'+html.escape(title)+'</h1><p>'+html.escape(notice)+'</p><p>'+html.escape(note)+'</p></header>'+''.join(blocks)+'<footer>每张图的源数据与显示编号对应表见 figures/source_data；绘图脚本、参数和输入SHA256见figure_manifest.json。参考 <a href="https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/">Nature图形指南</a>。图形文件自动生成后仍需最终人工版面检查。</footer></main></html>'
    (out/'图文结果报告.html').write_text(page,encoding='utf-8');(out/'图文结果报告.md').write_text('\n\n'.join(md)+'\n',encoding='utf-8')
    manifest=dict(status='synthetic_layout_preview' if a.synthetic_preview else 'generated_pending_visual_review',created_utc=datetime.now(timezone.utc).isoformat(),job_id=os.environ.get('SLURM_JOB_ID'),pilot=a.pilot,matplotlib=matplotlib.__version__,numpy=np.__version__,font=font,upstream_acceptance_sha256=sha(acceptance),input_sha256={str(v.relative_to(root)):sha(v) for v in [*files.values(),frozen,*acceptance_inputs]},script_sha256={p.name:sha(p) for p in [Path(__file__),HERE/'pilot_function_figures.py',HERE/'pilot10_catalog_acceptance.py']},figures=figures,outputs={str(p.relative_to(out)):sha(p) for p in [*figdir.rglob('*'),out/'图文结果报告.html',out/'图文结果报告.md'] if p.is_file()})
    (figdir/'figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(dict(status=manifest['status'],figures=len(figures),output=str(out/'图文结果报告.html')),ensure_ascii=False))
if __name__=='__main__':main()
