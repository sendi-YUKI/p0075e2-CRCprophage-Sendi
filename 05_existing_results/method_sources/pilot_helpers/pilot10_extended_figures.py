"""Two supplementary pilot figures; verified inputs, no new statistical tests.

Production must run in Slurm. Synthetic fixtures must be marked in both upstream
receipts and live below figure_style_preview. Existing five figures are untouched.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import html
import json
import math
import os
import sys

HERE = Path(__file__).resolve().parent
if (HERE / 'figure_preview_dependencies').is_dir():
    sys.path.insert(0, str(HERE / 'figure_preview_dependencies'))
COLORS = {'CRC': '#D55E00', 'control': '#0072B2', 'reference': '#009E73',
          'metagenome': '#8064A2', 'unknown': '#D9DEE3', 'ink': '#172B3A'}
SOURCE_LABELS = {'legacy_bacterial_genome': 'Reference genomes',
                 'metagenome_assembly': 'Stool metagenomes'}


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read_table(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def export(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row)) or ['status']
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    if value is None or str(value).strip().lower() in {'', 'na', 'nan', 'none', 'unknown', 'null'}:
        return None
    value = float(value)
    if not math.isfinite(value):
        return None
    return value


def state(value):
    text = str(value).strip().lower()
    if text in {'true', 'yes', '1'}:
        return True
    if text in {'false', 'no', '0'}:
        return False
    return None


def verify_receipt(path, expected_kind, synthetic):
    receipt = read_json(path)
    if receipt.get('status') != 'completed':
        raise ValueError(f'Upstream stage not accepted: {path}')
    kind_key = 'analysis_type' if expected_kind == 'prophage_evidence' else 'analysis_status'
    if receipt.get(kind_key) != expected_kind:
        raise ValueError(f'Unexpected upstream analysis identity: {path}')
    if bool(receipt.get('synthetic_preview', False)) != synthetic:
        raise ValueError('Synthetic and production inputs must never mix')
    outputs = receipt.get('outputs', {})
    if not outputs:
        raise ValueError(f'Missing upstream output hashes: {path}')
    for relative, digest in outputs.items():
        target = (path.parent / relative).resolve()
        if not target.is_relative_to(path.parent.resolve()):
            raise ValueError(f'Unsafe upstream output path: {relative}')
        if not target.is_file() or sha(target) != digest:
            raise ValueError(f'Upstream output hash mismatch: {target}')
    for source, digest in receipt.get('input_sha256', {}).items():
        source = Path(source)
        if not source.is_absolute():
            source = path.parent / source
        if not source.is_file() or sha(source) != digest:
            raise ValueError(f'Upstream input hash mismatch: {source}')
    return receipt


def panel(ax, letter, title):
    ax.set_title(title, loc='left', fontsize=8, pad=12)
    ax.text(-.12, 1.04, letter, transform=ax.transAxes, fontweight='bold', fontsize=11)


def empty(ax, text):
    ax.axis('off')
    ax.text(.5, .5, text, transform=ax.transAxes, ha='center', va='center', color='#526273')


def evidence_state(row, key):
    if key == 'predicted_provirus':
        kind = row.get('candidate_type', '')
        return None if not kind else kind == 'provirus_locus'
    if key == 'two_sided_flank_geometry':
        left, right = number(row.get('left_flank_bp')), number(row.get('right_flank_bp'))
        return None if left is None or right is None else left > 0 and right > 0
    return state(row.get(key))


def quality_label(value):
    value = str(value).strip().lower().replace('_', '-').replace(' ', '-')
    return {'complete': 'Complete', 'high-quality': 'High quality', 'high': 'High quality',
            'medium-quality': 'Medium quality', 'medium': 'Medium quality',
            'low-quality': 'Low quality', 'low': 'Low quality'}.get(value, 'Undetermined / NA')


def figure_evidence(c):
    plt, np, rows = c['plt'], c['np'], c['evidence']
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.3), layout='constrained',
                             gridspec_kw={'width_ratios': [1.25, 1]})
    dims = [('predicted_provirus', 'Predicted provirus locus'),
            ('two_sided_flank_geometry', 'Flanks on both sides (geometry)'),
            ('two_sided_host_marker_evidence', 'Host markers on both sides'),
            ('checkv_provirus', 'CheckV provirus signal')]
    quality = ['Complete', 'High quality', 'Medium quality', 'Low quality', 'Undetermined / NA']
    aggregate = []
    sources = list(SOURCE_LABELS)
    for axis, letter, title in zip(axes, 'ab', ['Separate evidence dimensions', 'CheckV quality categories']):
        panel(axis, letter, title)
    if not rows:
        for axis in axes:
            empty(axis, 'No occurrence evidence available')
    else:
        for x, source in enumerate(sources):
            subset = [row for row in rows if row['source_kind'] == source]
            color = COLORS['reference' if x == 0 else 'metagenome']
            for y, (key, label) in enumerate(dims):
                values = [evidence_state(row, key) for row in subset]
                positive, negative = values.count(True), values.count(False)
                unknown = values.count(None)
                known = positive + negative
                if known:
                    axes[0].scatter(x, y, s=25 + 250 * positive / known, c=color, alpha=.75, edgecolors='white', linewidths=.5)
                    text = f'{positive}/{known}'
                else:
                    axes[0].scatter(x, y, s=25, c='#8D99A3', marker='x', linewidths=.9)
                    text = 'NA'
                if unknown:
                    text += f'\n+{unknown} NA'
                axes[0].text(x + .17, y, text, va='center', fontsize=6)
                aggregate.append(dict(source_kind=source,dimension=key,total_occurrences=len(subset),positive=positive,negative=negative,unknown=unknown,assessed_denominator=known))
            counts = Counter(quality_label(row.get('checkv_quality')) for row in subset)
            for y, label in enumerate(quality):
                count = counts[label]
                if subset:
                    axes[1].scatter(x, y, s=18 + 250 * count / len(subset),
                                    facecolors=color if count else 'white', edgecolors=color, linewidths=.6, alpha=.8)
                    axes[1].text(x + .17, y, str(count), va='center', fontsize=6.5)
                else:
                    axes[1].scatter(x, y, s=25, c='#8D99A3', marker='x')
                    axes[1].text(x + .17, y, 'NA', va='center', fontsize=6.5)
                aggregate.append(dict(source_kind=source,dimension='checkv_quality:'+label,total_occurrences=len(subset),positive=count,negative='',unknown='',assessed_denominator=len(subset)))
        for ax, labels in zip(axes, [[item[1] for item in dims], quality]):
            ax.set_xlim(-.30, 1.62)
            ax.set_ylim(len(labels)-.5, -.6)
            ax.set_yticks(range(len(labels)), labels, fontsize=7)
            ax.set_xticks(range(2), [SOURCE_LABELS[s]+'\n('+str(sum(row['source_kind']==s for row in rows))+' occurrences)' for s in sources], fontsize=6.5)
            ax.tick_params(axis='both', length=0)
            for spine in ax.spines.values():
                spine.set_visible(False)
        axes[0].set_xlabel('Positive / assessed occurrences; NA shown separately', labelpad=12, fontsize=6.5)
        axes[1].set_xlabel('Occurrence counts; size = within-source fraction', labelpad=12, fontsize=6.5)
    c['save'](fig, 'Fig06_evidence_quality')
    c['export']('Fig06_occurrence_evidence', rows)
    c['export']('Fig06_aggregate_dot_matrix', aggregate)
    return dict(stem='Fig06_evidence_quality', title='按来源查看 prophage 证据与序列质量',
        caption='左图的每一行是独立证据维度，数字为阳性数／可评估数，另列 NA；不是从低到高的真伪等级，也不能将行数相加。两侧几何侧翼只表示序列空间，不能自动称为细菌侧翼。CheckV 未报告 provirus 不等于否定 geNomad 的预测。右图按原始出现实例统计质量类别；无法确定质量单独保留，不视为零质量。参考端与宏基因组端各有自己的分母，同一病毒可在多个来源出现，实例数不等于独立 vOTU 数。',
        selection='All accepted occurrence rows, grouped by source; evidence dimensions and quality shown separately.')


def selected_votus(rows):
    def key(row):
        total = (number(row.get('crc_detection_n')) or 0) + (number(row.get('control_detection_n')) or 0)
        found = (number(row.get('crc_detected')) or 0) + (number(row.get('control_detected')) or 0)
        return (-(found / total) if total > 0 else float('inf'), row['votu_id'])
    return sorted(rows, key=key)[:16]


def qtext(value):
    value = number(value)
    return 'NA' if value is None else ('<.001' if value < .001 else f'{value:.3f}' if value < .01 else f'{value:.2f}')


def figure_comparison(c):
    plt, np = c['plt'], c['np']
    selected, donors = c['selected'], c['donors']
    fig, axes = plt.subplots(1, 3, figsize=(10.6, max(4.0, 2.3 + len(selected)*.20)),
                             layout='constrained', gridspec_kw={'width_ratios': [1, 1, 2.15]})
    for ax, letter, title in zip(axes, 'abc', ['Detection difference', 'Abundance effect', 'Individual donor abundance']):
        panel(ax, letter, title)
    source = []
    if not selected:
        for ax in axes:
            empty(ax, 'No vOTUs available for display')
    else:
        for ax, effect, qcol, xlabel in [(axes[0], 'detection_fraction_difference','fisher_q_global','CRC - control detection fraction'),
                                        (axes[1], 'cliffs_delta','abundance_q_global',"Cliff's delta (CRC vs control)")]:
            ax.axvline(0, color='#8D99A3', linestyle='--', linewidth=.7, zorder=1)
            for i, row in enumerate(selected):
                value = number(row.get(effect))
                if value is not None:
                    ax.scatter(value, i, s=30, c=COLORS['CRC'] if value > 0 else COLORS['control'] if value < 0 else '#526273', edgecolors='white', linewidths=.5, zorder=3)
                else:
                    ax.text(0, i, 'NA', ha='center', va='center', color='#8D99A3', fontsize=6)
                ax.text(1.25, i, qtext(row.get(qcol)), ha='center', va='center', fontsize=6)
            ax.set_xlim(-1.08, 1.53)
            ax.set_xticks([-1, 0, 1])
            ax.set_ylim(len(selected)-.5, -.8)
            ax.set_yticks(range(len(selected)), ['V'+str(i+1).zfill(2) for i in range(len(selected))] if ax is axes[0] else [])
            ax.set_xlabel(xlabel, fontsize=6.5, labelpad=9)
            ax.text(1.25, -.75, 'Global q', ha='center', va='bottom', fontsize=6)
        lookup = {(r['votu_id'],r['biological_sample_id']):r for r in c['donor_values']}
        ax = axes[2]
        positive = []
        for i, stat in enumerate(selected):
            for j, donor in enumerate(donors):
                row = lookup.get((stat['votu_id'],donor))
                value = number(row.get('abundance_rpkm')) if row and state(row.get('abundance_valid')) is True else None
                if value is None:
                    ax.scatter(j, i, c='#A5AEB7', marker='x', s=17, linewidths=.75)
                elif value == 0:
                    ax.scatter(j, i, facecolors='white', edgecolors='#596B79', s=21, linewidths=.65)
                else:
                    positive.append((j,i,math.log10(1+value)))
                source.append(dict(row or {},votu_id=stat['votu_id'],biological_sample_id=donor,display_alias='V'+str(i+1).zfill(2),display_log10_one_plus_rpkm='' if value is None else math.log10(1+value),display_state='NA' if value is None else 'zero' if value==0 else 'positive'))
        if positive:
            im = ax.scatter([p[0] for p in positive],[p[1] for p in positive],c=[p[2] for p in positive],s=25,cmap='viridis',vmin=0)
            fig.colorbar(im,ax=ax,orientation='horizontal',pad=.025,shrink=.75,aspect=28,label='Positive values: log10(1 + fragment RPKM)')
        for j, donor in enumerate(donors):
            ax.scatter(j,-.98,c=COLORS[c['groups'][donor]],s=22,marker='s',clip_on=False)
        ax.set_xticks(range(len(donors)),donors,rotation=90,fontsize=5.8)
        ax.set_yticks([])
        ax.set_xlim(-.6,max(.6,len(donors)-.4))
        ax.set_ylim(len(selected)-.5,-.8)
        ax.spines['left'].set_visible(False)
        from matplotlib.lines import Line2D
        handles = [Line2D([],[],marker='s',ls='',markersize=4,color=COLORS[g],label=g) for g in ['CRC','control']]
        handles += [Line2D([],[],marker='o',ls='',markersize=4,markerfacecolor='white',color='#596B79',label='Measured zero'),Line2D([],[],marker='x',ls='',markersize=4,color='#A5AEB7',label='Not assessable')]
        fig.legend(handles=handles,loc='outside lower center',ncol=4,fontsize=6)
    c['save'](fig,'Fig07_exploratory_comparison')
    c['export']('Fig07_selected_votu_statistics',[dict(row,display_alias='V'+str(i+1).zfill(2)) for i,row in enumerate(selected)])
    c['export']('Fig07_individual_abundance',source)
    return dict(stem='Fig07_exploratory_comparison',title='检出与丰度的探索性组间比较',
        caption='最多展示 16 个 vOTU，按全部供者合并的检出比例降序、再按 vOTU 编号排序，未按 P 值、q 值或 CRC 效应选择。左侧为 CRC 减对照的检出比例差，中间为丰度的 Cliff’s delta；没有置信区间输入，因此不画误差线。两列 q 值均来自全部 vOTU 的检出与丰度两个端点合并后的全局 BH 校正，不采用逐端点 q 值。右侧每列一名供者，空心圆是有效测得的零，灰叉是不可评估，彩点是正丰度。分母和缺失原因见源表。本批发现与测量使用同一批样本，小样本差异只能用于提出假说；这不是独立验证、宿主校正效应或患者功能基因表达。',
        selection='Pooled detection prevalence descending, then vOTU ID; at most 16, never ranked by P/q/effect.')


def write_report(c, figures):
    out, figdir = c['out'], c['figdir']
    title='CRC-PHIRE · 先导补充图文报告'
    notice='合成排版预览，不是研究结果。' if c['synthetic'] else '真实数据作图；探索性比较不代表正式疾病结论。'
    ntest=sum(number(r.get('fisher_q_global')) is not None for r in c['stats'])+sum(number(r.get('abundance_q_global')) is not None for r in c['stats'])
    nlow=sum(number(r.get(col)) is not None and number(r.get(col))<.05 for r in c['stats'] for col in ['fisher_q_global','abundance_q_global'])
    summary=f"证据表包含 {len(c['evidence'])} 个来源实例；统计表包含 {len(c['stats'])} 个 vOTU。可报告全局 q 值的端点共 {ntest} 个，其中 q＜0.05 的端点为 {nlow} 个；同一 vOTU 可以贡献两个端点，不能把这个数当作独立候选数量。"
    links=[('原五图图文报告','图文结果报告.html'),('全部探索性比较','exploratory/votu_differential.tsv'),('逐人测量与缺失原因','exploratory/donor_analysis_values.tsv'),('完整来源证据表','prophage_evidence/occurrence_evidence.tsv')]
    if (out/'exploratory/differential_function_links.tsv').exists():
        links.append(('差异结果与功能注释连接','exploratory/differential_function_links.tsv'))
    md=[f'# {title}',notice,summary,'本页补充原有五图，不覆盖其文件。证据维度用于解释序列来源和结构；统计差异用于探索，二者不互相替代。',' · '.join(f'[{label}]({url})' for label,url in links)]
    blocks=[]
    for i,f in enumerate(figures,6):
        stem=f['stem']
        md += [f'## Figure {i}：{f["title"]}',f'![{f["title"]}](extended_figures/{stem}.png)',f['caption'],f'展示规则：{f["selection"]}',f'[SVG](extended_figures/{stem}.svg) · [PDF](extended_figures/{stem}.pdf)']
        blocks.append('<section><div class="eyebrow">FIGURE '+str(i).zfill(2)+'</div><h2>'+html.escape(f['title'])+'</h2><img src="extended_figures/'+stem+'.svg" alt="'+html.escape(f['title'])+'"><p>'+html.escape(f['caption'])+'</p><p class="rule">'+html.escape(f['selection'])+'</p><p>'+ ' · '.join(f'<a href="extended_figures/{stem}.{ext}">{ext.upper()}</a>' for ext in ['svg','pdf','png'])+'</p></section>')
    if (out/'exploratory/differential_function_links.tsv').exists():
        md += ['功能连接表把序列本身的 PHROGs／系统／KO 注释连接到 vOTU 探索性统计；它不是功能差异检验，不表示患者功能基因丰度、携带或表达。']
    aliases=[dict(display_alias='V'+str(i+1).zfill(2),votu_id=r['votu_id'],reference_id=r.get('reference_id','')) for i,r in enumerate(c['selected'])]
    c['export']('votu_display_aliases',aliases,fields=['display_alias','votu_id','reference_id'])
    alias_html='<details><summary>显示编号与原始 vOTU 对应</summary><table><tr><th>图中编号</th><th>vOTU</th><th>参考编号</th></tr>'+''.join('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in ['display_alias','votu_id','reference_id'])+'</tr>' for r in aliases)+'</table></details>'
    style='body{margin:0;background:#f4f6f7;color:#192d3b;font:16px/1.8 system-ui,"Microsoft YaHei",sans-serif}main{max-width:1180px;margin:auto;padding:36px 24px}header{padding:18px 0 28px;border-bottom:2px solid #192d3b}h1{font-size:30px;line-height:1.35}h2{font-size:23px}.eyebrow{font-size:12px;letter-spacing:2px;color:#526575}section{background:white;padding:28px;margin:26px 0;border:1px solid #dce3e7;border-radius:8px}img{display:block;width:100%;height:auto;margin:22px auto}.rule,footer{font-size:13px;color:#526575}a{color:#0072b2}details{overflow:auto}table{border-collapse:collapse;font-size:12px}td,th{padding:7px;border-bottom:1px solid #dce3e7;text-align:left;overflow-wrap:anywhere}@media(max-width:650px){main{padding:18px 12px}section{padding:16px}h1{font-size:25px}}'
    page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+title+'</title><style>'+style+'</style><main><header><div class="eyebrow">CRC-PHIRE / SUPPLEMENTARY PILOT FIGURES</div><h1>'+title+'</h1><p>'+notice+'</p><p>'+summary+'</p><p>'+' · '.join('<a href="'+url+'">'+label+'</a>' for label,url in links)+'</p></header>'+''.join(blocks)+alias_html+'<footer><p>全部显示数据见 extended_figures/source_data；输入与输出校验值见 figure_manifest.json。功能目录或功能组合的探索结果，不等于患者中的功能基因表达或机制实验证据。</p></footer></main></html>'
    (out/'先导补充图文报告.md').write_text('\n\n'.join(md)+'\n',encoding='utf-8')
    (out/'先导补充图文报告.html').write_text(page,encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=HERE.parent)
    parser.add_argument('--synthetic-preview',action='store_true')
    args=parser.parse_args()
    root=args.root.resolve(); out=root/'12_reports/pilot10'
    if args.synthetic_preview:
        if 'figure_style_preview' not in root.parts:
            raise ValueError('Synthetic previews require a separate figure_style_preview root')
    elif not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('Production figures must run through Slurm')
    evidence_manifest=out/'prophage_evidence/evidence_manifest.json'
    stats_manifest=out/'exploratory/exploratory_manifest.json'
    er=verify_receipt(evidence_manifest,'prophage_evidence',args.synthetic_preview)
    sr=verify_receipt(stats_manifest,'completed_exploratory_5_vs_5',args.synthetic_preview)
    files={'evidence':out/'prophage_evidence/occurrence_evidence.tsv',
           'stats':out/'exploratory/votu_differential.tsv',
           'donor_values':out/'exploratory/donor_analysis_values.tsv'}
    for name,path in files.items():
        receipt=er if name=='evidence' else sr
        if path.name not in receipt['outputs']:
            raise ValueError(f'Unverified figure input: {path}')
    data={key:read_table(path) for key,path in files.items()}
    if len({r['votu_id'] for r in data['stats']}) != len(data['stats']):
        raise ValueError('Duplicated vOTU statistics')
    for row in data['stats']:
        for col in ['fisher_q_global','abundance_q_global']:
            value=number(row.get(col))
            if value is not None and not 0 <= value <= 1:
                raise ValueError(f'Invalid global q value: {row["votu_id"]}/{col}')
    if len({(r['votu_id'],r['biological_sample_id']) for r in data['donor_values']}) != len(data['donor_values']):
        raise ValueError('Duplicated donor x vOTU measurement')
    if any(row['source_kind'] not in SOURCE_LABELS for row in data['evidence']):
        raise ValueError('Unexpected source kind; do not silently omit it')
    groups={}
    for row in data['donor_values']:
        donor, group=row['biological_sample_id'],row['group']
        if group not in {'CRC','control'} or donor in groups and groups[donor]!=group:
            raise ValueError('Invalid or conflicting donor groups')
        groups[donor]=group
        value=number(row.get('abundance_rpkm'))
        if value is not None and value<0:
            raise ValueError('Negative RPKM')
    figdir=out/'extended_figures'
    if (figdir/'figure_manifest.json').exists():
        raise ValueError('Preserve accepted figures; use a new reviewed output version to redraw')
    if not args.synthetic_preview and not (out/'图文结果报告.html').is_file():
        raise ValueError('Original five-figure report must exist; no broken report link')
    figdir.mkdir(parents=True,exist_ok=True); (figdir/'source_data').mkdir(exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR',str(root/'92_tmp/extended_figure_mplconfig'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import font_manager
    font='Arial' if any(f.name=='Arial' for f in font_manager.fontManager.ttflist) else 'DejaVu Sans'
    plt.rcParams.update({'font.family':font,'font.size':7,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.65,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white','text.color':COLORS['ink'],'axes.labelcolor':COLORS['ink'],'savefig.dpi':300})
    def save(fig,stem):
        if args.synthetic_preview:
            fig.suptitle('SYNTHETIC LAYOUT TEST - NOT STUDY RESULTS',fontsize=8,color='#AB3900')
        for ext in ['svg','pdf','png']:
            fig.savefig(figdir/(stem+'.'+ext),bbox_inches='tight',pad_inches=.12)
        plt.close(fig)
    c=dict(data,root=root,out=out,figdir=figdir,plt=plt,np=np,synthetic=args.synthetic_preview,groups=groups,
           donors=sorted(groups,key=lambda s:(groups[s]!='CRC',s)),selected=selected_votus(data['stats']),
           save=save,export=lambda stem,rows,fields=None:export(figdir/'source_data'/(stem+'.tsv'),rows,fields))
    figures=[figure_evidence(c),figure_comparison(c)]
    write_report(c,figures)
    all_inputs=[evidence_manifest,stats_manifest,*files.values()]
    manifest=dict(status='synthetic_layout_preview' if args.synthetic_preview else 'generated_pending_visual_review',
        synthetic_preview=args.synthetic_preview,created_utc=datetime.now(timezone.utc).isoformat(),job_id=os.environ.get('SLURM_JOB_ID'),
        matplotlib=matplotlib.__version__,numpy=np.__version__,font=font,
        input_sha256={str(p.relative_to(root)):sha(p) for p in all_inputs},script_sha256={Path(__file__).name:sha(Path(__file__))},
        figures=figures,statistical_rule='Only global BH q values; no tests recomputed; no confidence intervals manufactured',
        output_sha256={str(p.relative_to(out)):sha(p) for p in [*figdir.rglob('*'),out/'先导补充图文报告.md',out/'先导补充图文报告.html'] if p.is_file()})
    (figdir/'figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(dict(status=manifest['status'],figures=len(figures),output=str(out/'先导补充图文报告.html')),ensure_ascii=False))


if __name__=='__main__':
    main()
