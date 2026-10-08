import pathlib,csv,json,shutil,hashlib,argparse,collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
R=pathlib.Path('/srv/CRC-PHIRE/analysis1');A=R/'14_postpilot_adjustment';O=A/'07_figures';D=O/'source_data'
O.mkdir(exist_ok=True);D.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
COL={'CRC':'#D97739','control':'#347EA5'}
def table(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def savefig(fig,name):
    for ext in ['png','pdf','svg']:fig.savefig(O/(name+'.'+ext),dpi=220,bbox_inches='tight')
    plt.close(fig)
def source(p):shutil.copyfile(p,D/p.name)
def panel(ax,letter,title):ax.set_title(letter+'  '+title,loc='left',fontsize=11,fontweight='bold',pad=12)
def core():
    audit=A/'04_measurement_audit';rep=A/'05_replication_activity';cal=audit/'technical_calibration'
    for path in [audit/'ambiguity_fragment_classes.tsv',audit/'threshold_summary.tsv',cal/'calibration_results.tsv',cal/'selected_technical_examples.tsv',rep/'locus_diagnostics.tsv',rep/'coverage_bins.tsv']:source(path)
    ar=table(audit/'ambiguity_fragment_classes.tsv');tot=collections.Counter()
    for r in ar:tot[r['ambiguity_class']]+=int(r['unique_within_readset_fragments'])
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    keys=['host_only','multiple_targets','target_host','same_target_locations_or_saturation'];labels=['Host-background only','Between viral targets','Viral target + host','Within-target / saturated']
    aa=ax[0,0];aa.scatter([tot[k] for k in keys],range(4),s=75,c=['#9CA3AF','#7566A7','#347EA5','#65A995'])
    for i,k in enumerate(keys):aa.text(tot[k]*1.12,i,f'{tot[k]:,}',va='center',fontsize=8)
    aa.set_yticks(range(4),labels);aa.set_xscale('log');aa.set_xlim(500,max(tot.values())*4);aa.invert_yaxis();aa.set_xlabel('Ambiguous fragments (counted once per readset)');panel(aa,'A','Where ambiguity occurs')
    tr=table(audit/'threshold_summary.tsv');states=['detected','not_detected','not_detected_at_resolved_signal','not_resolvable'];matrix=np.array([[next(int(r['sample_votu_cells']) for r in tr if float(r['threshold'])==t and r['group']==g and r['status']==s) for s in states] for t,g in [(.75,'CRC'),(.75,'control'),(.30,'CRC'),(.30,'control')]])
    aa=ax[0,1];im=aa.imshow(matrix,cmap='Blues',aspect='auto')
    for i in range(4):
        for j in range(4):aa.text(j,i,str(matrix[i,j]),ha='center',va='center',color='white' if matrix[i,j]>1000 else '#243447')
    aa.set_xticks(range(4),['Detected','Below rule','Mixed unresolved','Ambiguous only'],rotation=25,ha='right');aa.set_yticks(range(4),['CRC / 0.75','Control / 0.75','CRC / 0.30','Control / 0.30']);panel(aa,'B','Sample x vOTU states (2,365 cells per row)')
    sim=table(cal/'calibration_results.tsv');targets=sorted({r['target'] for r in sim});alias={t:'T'+str(i+1).zfill(2) for i,t in enumerate(targets)}
    aa=ax[1,0]
    for layout,offset,color in [('SE',-.1,'#347EA5'),('PE',.1,'#D97739')]:
        y=[float(next(r['breadth'] for r in sim if r['target']==t and r['scenario']=='positive' and r['layout']==layout)) for t in targets]
        aa.scatter(np.arange(len(targets))+offset,y,s=36,c=color,label=layout+' known positive')
    aa.axhline(.75,c='#555',ls='--',lw=.8);aa.axhline(.3,c='#888',ls=':',lw=.8);aa.set_xticks(range(len(targets)),[alias[t] for t in targets],rotation=45);aa.set_ylim(-.03,1.06);aa.set_ylabel('Uniquely assigned coverage breadth');aa.legend(frameon=False,fontsize=8,loc='lower right');panel(aa,'C','Bounded technical challenge (synthetic)')
    aa=ax[1,1]
    for layout,marker in [('SE','o'),('PE','s')]:
        vals=[float(next(r['breadth'] for r in sim if r['target']==t and r['scenario']=='negative' and r['layout']==layout)) for t in targets]
        aa.scatter(range(len(targets)),vals,s=45,marker=marker,facecolors='none',edgecolors='#347EA5' if layout=='SE' else '#D97739',label=layout+' selected target absent')
    aa.set_xticks(range(len(targets)),[alias[t] for t in targets],rotation=45);aa.set_ylim(-.03,1.03);aa.axhline(.3,c='#888',ls=':',lw=.8);aa.set_ylabel('Coverage wrongly attributed to selected target');aa.legend(frameon=False,fontsize=8);panel(aa,'D','Near-neighbour and host-background challenge')
    fig.suptitle('Measurement diagnosis: ambiguity, breadth and known-source checks',fontsize=15,fontweight='bold');savefig(fig,'Fig02_measurement_diagnosis')
    (D/'technical_target_aliases.json').write_text(json.dumps(alias,indent=2))
    lr=table(rep/'locus_diagnostics.tsv');bins=table(rep/'coverage_bins.tsv');active=[r for r in lr if r['replication_call']=='replication_compatible'];qualified=[r for r in lr if r['membership_status'] in ['member','representative'] and r['replication_call']!='replication_compatible']
    near=sorted(qualified,key=lambda r:(abs(float(r['tool_prophage-host_ratio'])-1),r['source_occurrence_id']));selection=active+next(([r] for r in near if r['group']=='CRC'),[])+next(([r] for r in near if r['group']=='control'),[])
    fig,axes=plt.subplots(4,1,figsize=(12,10),layout='constrained')
    for i,(aa,r) in enumerate(zip(axes,selection)):
        rows=[b for b in bins if b['source_occurrence_id']==r['source_occurrence_id']];x=np.array([(int(b['start0'])+int(b['end0']))/2000 for b in rows]);y=np.array([float(b['native_depth']) if b['native_depth']!='NA' else np.nan for b in rows]);z=np.array([float(b['strict_depth']) if b['strict_depth']!='NA' else np.nan for b in rows]);a,b=int(r['start0'])/1000,int(r['end0'])/1000
        aa.axvspan(a,b,color='#EEE7D8',alpha=.9,label='Predicted provirus');aa.plot(x,y,c='#347EA5',lw=1,label='Native span coverage');aa.plot(x,z,c='#D97739',lw=.8,alpha=.8,label='Strict CIGAR-block coverage');aa.set_ylabel('Depth (x)');aa.set_xlabel('Original contig position (kb)');aa.set_xlim(x.min(),x.max());aa.set_ylim(bottom=0)
        label='Replication candidate' if r['replication_call']=='replication_compatible' else 'Near-unity example'
        panel(aa,chr(65+i),f"{r['biological_sample_id']} ({r['group']}) | {label} | ratio={float(r['tool_prophage-host_ratio']):.2f}, local={float(r['local_ratio']):.2f}")
        aa.text(.99,.9,f"Provirus {(b-a):.1f} kb | {r['membership_status'].replace('excluded_viral_eligibility_ineligible','outside original qualified catalog')}",ha='right',transform=aa.transAxes,fontsize=8)
        if i==0:aa.legend(frameon=False,ncol=3,loc='upper left',fontsize=8)
    fig.suptitle('Original assembly context supports specific replication candidates',fontsize=15,fontweight='bold');savefig(fig,'Fig03_locus_coverage')
    (D/'coverage_example_selection.json').write_text(json.dumps({'selection':'both tool-positive loci; closest-to-one qualified locus per group; no disease P selection','occurrence_ids':[r['source_occurrence_id'] for r in selection]},indent=2))
    fig,axes=plt.subplots(1,3,figsize=(14,4.8),layout='constrained');aa=axes[0]
    for group in ['CRC','control']:
        rr=[r for r in lr if r['group']==group];aa.scatter([float(r['tool_host_mean_cov']) for r in rr],[float(r['tool_prophage_mean_cov']) for r in rr],s=[65 if r['replication_call']=='replication_compatible' else 16 for r in rr],c=COL[group],alpha=.65,label=group)
    lim=max(max(float(r['tool_host_mean_cov']),float(r['tool_prophage_mean_cov'])) for r in lr)*1.3;aa.plot([1,lim],[1,lim],c='#555',lw=.8);aa.plot([1,lim/2],[2,lim],c='#999',ls='--',lw=.8);aa.set_xscale('log');aa.set_yscale('log');aa.set_xlabel('Host-region mean depth');aa.set_ylabel('Provirus mean depth');aa.legend(frameon=False);panel(aa,'A','382 coverage-evaluable loci')
    aa=axes[1];samples=sorted({r['biological_sample_id'] for r in lr},key=lambda sid:(next(r['group'] for r in lr if r['biological_sample_id']==sid)!='CRC',sid))
    for i,sid in enumerate(samples):
        rr=[r for r in lr if r['biological_sample_id']==sid];ratios=[float(r['tool_prophage-host_ratio']) for r in rr];jitter=np.linspace(-.15,.15,len(rr));aa.scatter(i+jitter,ratios,c=COL[rr[0]['group']],s=12,alpha=.65);aa.text(i,3.75,str(len(rr)),ha='center',fontsize=8)
    aa.axhline(1,c='#777',lw=.7);aa.axhline(2,c='#777',ls='--',lw=.7);aa.set_xticks(range(len(samples)),samples,rotation=65,ha='right',fontsize=7);aa.set_ylabel('Provirus / host coverage ratio');aa.set_ylim(0,4);panel(aa,'B','Within-donor distribution (n above)')
    aa=axes[2]
    for typ,color,marker in [('strict_ratio','#347EA5','o'),('local_ratio','#D97739','x')]:aa.scatter([float(r['tool_prophage-host_ratio']) for r in lr],[float(r[typ]) for r in lr],s=15,c=color,alpha=.55,marker=marker,label='Strict mapping' if typ=='strict_ratio' else 'Local host denominator')
    mx=max(float(r['local_ratio']) for r in lr)*1.1;aa.plot([0,mx],[0,mx],c='#888',lw=.8);aa.set_xlim(0,mx);aa.set_ylim(0,mx);aa.set_xlabel('Primary coverage ratio');aa.set_ylabel('Sensitivity coverage ratio');aa.legend(frameon=False,fontsize=8);panel(aa,'C','Finite method sensitivity')
    fig.suptitle('Replication branch: coverage feasibility and candidate stability',fontsize=15,fontweight='bold');savefig(fig,'Fig04_replication_feasibility')
def host():
    assert (A/'03_host_profile/completed.json').is_file()
    rows=table(A/'03_host_profile/donor_species.tsv');identity=table(A/'01_manifests/donor_input_identity.tsv');groups={r['sample']:r['group'] for r in identity};samples=sorted(groups,key=lambda x:(groups[x]!='CRC',x));taxa=collections.Counter();ab={}
    for r in rows:
        if r['taxon_level']=='species' and r['clade_name'].startswith('k__Bacteria'):
            t=r['clade_name'].split('|')[-1][3:];v=float(r['relative_abundance_percent']);taxa[t]+=v;ab[(r['biological_sample_id'],t)]=v
    priority=['Bacteroides_fragilis','Bacteroides_hominis','Fusobacterium_animalis','Fusobacterium_nucleatum','Parvimonas_micra','Peptostreptococcus_anaerobius']
    selected=[t for t in priority if t in taxa]+[t for t,_ in taxa.most_common(20) if t not in priority];selected=selected[:22]
    positive_min=min(v for (sid,t),v in ab.items() if t in selected and v>0)
    display_floor=10**(np.floor(np.log10(positive_min))-1)
    (D/'host_display_settings.json').write_text(json.dumps({'nonreported_display_floor_percent':display_floor,'floor_is_not_abundance':True,'positive_values_not_clipped':True},indent=2))
    fig,axes=plt.subplots(1,2,figsize=(13,7),gridspec_kw={'width_ratios':[1.4,1]},layout='constrained');aa=axes[0];matrix=np.array([[ab.get((sid,t),0) for sid in samples] for t in selected]);masked=np.ma.masked_less_equal(matrix,0);cm=plt.colormaps['viridis'].copy();cm.set_bad('#EEEEEE');im=aa.imshow(masked,norm=LogNorm(vmin=min(.001,positive_min),vmax=max(10,matrix.max())),cmap=cm,aspect='auto');aa.set_yticks(range(len(selected)),[t.replace('_',' ') for t in selected],fontsize=8);aa.set_xticks(range(10),samples,rotation=65,ha='right',fontsize=8)
    for tick,sid in zip(aa.get_xticklabels(),samples):tick.set_color(COL[groups[sid]])
    fig.colorbar(im,ax=aa,label='MetaPhlAn relative abundance (%)',shrink=.6);panel(aa,'A','Bacterial host context')
    aa=axes[1];focus=[t for t in priority if t in taxa];x=np.arange(len(focus))
    for group,offset in [('CRC',-.15),('control',.15)]:
        ids=[s for s in samples if groups[s]==group]
        for j,sid in enumerate(ids):aa.scatter(x+offset+(j-2)*.025,[(ab.get((sid,t),0) or display_floor) for t in focus],c=COL[group],s=26,alpha=.85,label=group if j==0 else None)
    aa.set_yscale('log');aa.set_xticks(x,[t.replace('_',' ') for t in focus],rotation=45,ha='right');aa.set_ylabel('Relative abundance (%)');aa.axhline(display_floor,c='#999',ls=':',lw=.7);aa.text(.02,.18,'Bottom row: not reported by profiler\nSpecies signal is not strain-level carriage',transform=aa.transAxes,fontsize=8);aa.legend(frameon=False);panel(aa,'B','Priority taxa: every donor shown')
    fig.suptitle('Host measurement complements viral sequence evidence',fontsize=15,fontweight='bold');savefig(fig,'Fig01_host_context');source(A/'03_host_profile/donor_species.tsv')
    with open(D/'host_heatmap_values.tsv','w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['sample','group','taxon','relative_abundance_percent'])
        for sid in samples:
            for t in selected:w.writerow([sid,groups[sid],t,ab.get((sid,t),0)])
p=argparse.ArgumentParser();p.add_argument('--part',choices=['core','host'],required=True);args=p.parse_args();globals()[args.part]()
manifest={'status':'generated_pending_visual_review','part':args.part,'matplotlib':matplotlib.__version__,'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.glob('Fig*') if p.is_file()},'synthetic_only_in_Fig02_panels_C_D':True}
(O/('manifest_'+args.part+'.json')).write_text(json.dumps(manifest,indent=2));print(json.dumps({'part':args.part,'figures':len(manifest['files'])}))
