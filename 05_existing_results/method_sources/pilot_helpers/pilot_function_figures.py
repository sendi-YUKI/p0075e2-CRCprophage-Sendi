"""Sequence-specific functional figures for completed CRC-PHIRE pilot reports.

Only observed annotation evidence is drawn. These figures make no host assignment,
AMG validation, gene synteny, disease association, or biological absence claim.
The caller owns input validation, output hashes and the Matplotlib backend.
"""
from collections import defaultdict
import re


INK = '#203746'
MUTED = '#6C7D87'
GRID = '#E2E8EC'
COLORS = ['#277F85', '#507BA8', '#946B98', '#CD8742']


def _true(value):
    return str(value).lower() in {'true', '1'}


def _activity(value):
    value = str(value).lower().replace('-', '').replace('_', '').replace(' ', '')
    if 'antidef' in value:
        return 'Anti-defense'
    if 'def' in value:
        return 'Defense'
    return None


def _sid(row):
    return row.get('viral_sequence_id') or row.get('contig_id') or ''


def _gene_ids(value):
    return [s for s in re.split(r'[;,|\s]+', str(value or '')) if s]


def _empty_panel(ctx, stem, title, detail, caption):
    fig, ax = ctx['plt'].subplots(figsize=(7.08, 2.2))
    ax.set_axis_off()
    ax.text(0, .90, title, transform=ax.transAxes, fontsize=11,
            fontweight='bold', color=INK)
    ax.text(0, .56, detail, transform=ax.transAxes, fontsize=9, color=MUTED,
            va='top', wrap=True)
    fig.subplots_adjust(left=.06, right=.97, top=.94, bottom=.08)
    ctx['save'](fig, stem)
    ctx['export'](stem, [{'status': 'no_displayable_evidence', 'reason': detail}])
    return {'stem': stem, 'title': title, 'caption': caption,
            'selection': 'No data were invented to fill an empty evidence panel.'}


def _evidence(ctx):
    tables = ctx['tables']
    genes = {r['gene_id']: r for r in tables.get('genes', [])}
    byseq = defaultdict(lambda: defaultdict(set))
    annotations = defaultdict(lambda: defaultdict(set))
    # function_hits contains rejected domains as well as accepted best hits.
    # The compact intersection table alone cannot establish acceptance.
    for row in tables.get('phrog_hits', []):
        if not (_true(row.get('accepted')) and _true(row.get('is_best'))):
            continue
        gid = row.get('gene_id', '')
        gene = genes.get(gid)
        feature = row.get('phrog', '')
        if gene and feature:
            byseq[_sid(gene)]['PHROG'].add(str(feature))
            annotations[gid]['PHROG'].add(str(feature))
    # Systems are called on each viral sequence's own ordered gene set.
    for row in tables.get('systems', []):
        sid = _sid(row)
        kind = _activity(row.get('activity', ''))
        feature = row.get('subtype') or row.get('type')
        if sid and kind and feature:
            byseq[sid][kind].add(str(feature))
            for gid in _gene_ids(row.get('all_gene_ids') or row.get('required_gene_ids')):
                if gid in genes and _sid(genes[gid]) == sid:
                    annotations[gid][kind].add(str(feature))
    for row in tables.get('features', []):
        sid = row.get('viral_sequence_id', '')
        if row.get('type') != 'KO_homology' or not sid or not row.get('feature'):
            continue
        feature = str(row['feature'])
        byseq[sid]['KO'].add(feature)
        for gid in _gene_ids(row.get('gene_id')):
            if gid in genes and _sid(genes[gid]) == sid:
                annotations[gid]['KO'].add(feature)
    return genes, byseq, annotations


def _dot_matrix(ctx, selected, evidence):
    stem = 'Fig04_sequence_function_evidence'
    title = 'Sequence-specific functional evidence'
    caption = ('每行是一条宏基因组来源病毒序列；圆点面积随该序列上不同功能标签数增加，'
               '数字给出实际数值。PHROG仅计accepted=true且is_best=true的同源家族；'
               '防御和反防御计不同系统亚型，KO计不同KO编号，四列的统计单位并不相同。'
               '“—”表示没有报告该类命中，“?”表示未评估或来源表缺失；均不能解释为生物学缺失。'
               '这些是序列自身的计算证据，KO不是已确认AMG，系统命中不是机制验证。'
               '不从相似参考序列转移注释，不进行组间功能显著性检验。')
    if not selected:
        return _empty_panel(ctx, stem, title,
                            'No broadly matched metagenomic sequence is available for display.', caption)
    plt = ctx['plt']
    fig, ax = plt.subplots(figsize=(7.08, max(3.1, .28 * len(selected) + 1.65)))
    cats = ['PHROG', 'Defense', 'Anti-defense', 'KO']
    labels = ['PHROG families', 'Defense subtypes', 'Anti-defense subtypes', 'KO identifiers']
    source_present = {
        'PHROG': 'phrog_hits' in ctx['tables'],
        'Defense': 'systems' in ctx['tables'],
        'Anti-defense': 'systems' in ctx['tables'],
        'KO': 'features' in ctx['tables'],
    }
    statuses = {_sid(r): r for r in ctx['tables'].get('function_status', [])}
    exported = []
    for y, sid in enumerate(selected):
        if y % 2 == 0:
            ax.axhspan(y - .5, y + .5, color='#F6F8FA', zorder=0)
        for x, cat in enumerate(cats):
            features = sorted(evidence.get(sid, {}).get(cat, set()))
            n = len(features)
            unassessed = cat != 'KO' and statuses.get(sid, {}).get('assessability') != 'assessable_observed_sequence'
            if not source_present[cat] or (unassessed and not n):
                ax.text(x, y, '?', ha='center', va='center', color=MUTED, fontsize=9)
                state = 'source_table_missing' if not source_present[cat] else 'sequence_function_not_assessable'
            elif n:
                area = min(460, 36 + 52 * n ** .5)
                ax.scatter(x, y, s=area, color=COLORS[x], alpha=.87,
                           edgecolors='white', linewidths=.7, zorder=3)
                ax.text(x, y, str(n), ha='center', va='center', fontsize=7.5,
                        color='white', fontweight='bold', zorder=4)
                state = 'observed_annotation_evidence'
            else:
                ax.text(x, y, '\u2014', ha='center', va='center', color='#9CA9B1', fontsize=8)
                state = 'no_reported_hit_not_biological_absence'
            exported.append(dict(viral_sequence_id=sid, display_alias=ctx['aliases'].get(sid, sid),
                                 category=cat, distinct_feature_count=n if source_present[cat] else '',
                                 feature_ids=';'.join(features), evidence_state=state,
                                 sequence_function_status=statuses.get(sid, {}).get('status', 'unknown'),
                                 unit='distinct_PHROG_family_or_system_subtype_or_KO_identifier'))
    ax.set_yticks(range(len(selected)), [ctx['aliases'].get(sid, sid) for sid in selected])
    ax.set_xticks(range(4), labels)
    ax.tick_params(axis='both', length=0, labelsize=8, colors=INK, pad=8)
    ax.xaxis.tick_top()
    ax.set_xlim(-.55, 3.55)
    ax.set_ylim(len(selected) - .5, -.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.suptitle(title, x=.08, y=.99, ha='left', fontsize=11, color=INK, fontweight='bold')
    fig.text(.08, .026, 'Own-sequence annotation evidence  |  \u2014 No reported hit  |  ? Not assessed / missing source',
             fontsize=7.5, color=MUTED)
    fig.subplots_adjust(left=.13, right=.97, top=1 - .72 / fig.get_figheight(), bottom=.12)
    ctx['save'](fig, stem)
    ctx['export'](stem, exported)
    return dict(stem=stem, title=title, caption=caption,
                selection='At most 16 broad-matched M sequences selected by pooled donor count, then length and ID; no CRC contrast ranking.')


def _gene_color(parts):
    if 'Anti-defense' in parts:
        return COLORS[2], 'Anti-defense system gene'
    if 'Defense' in parts:
        return COLORS[1], 'Defense system gene'
    if 'KO' in parts:
        return COLORS[3], 'KO homology'
    if 'PHROG' in parts:
        return COLORS[0], 'Accepted PHROG homology'
    return '#BCC7CF', 'No highlighted annotation'


def _gene_label(parts):
    for kind in ['Anti-defense', 'Defense', 'KO', 'PHROG']:
        if parts.get(kind):
            name = sorted(parts[kind])[0]
            name = str(name)
            if kind == 'PHROG':
                return 'PHROG ' + name.removeprefix('phrog_')
            return name if len(name) <= 25 else name[:22] + '...'
    return ''


def _gene_tracks(ctx, selected, genes, annotations):
    stem = 'Fig05_candidate_gene_context'
    title = 'Illustrative sequence-level gene context'
    caption = ('展示最多两条广泛匹配宏基因组病毒序列的自身基因坐标和方向；'
               '在先前按全体供者出现次数筛选的序列中按长度、编号排序选取，不按CRC与对照差异挑选。'
               '每个箭头是该序列实际预测的一个基因。颜色仅表示该基因自身的同源命中或系统成员证据，'
               '每条轨道最多标注六个基因；未着色不等于没有功能。'
               '同一基因可能有多类证据，所有类别保存在对应source-data表；颜色优先级为反防御、'
               '防御、KO、PHROG，仅用于显示。仅从系统表提供的实际gene ID连接系统成员。'
               '这些是未经过人工机制审查的展示实例，不是已确认实验候选；图中没有同源基因连线、'
               '序列共线性或宿主侧翼，不据此认定AMG、prophage整合或患者体内宿主。')
    seq = {r['viral_sequence_id']: r for r in ctx['tables'].get('seq', [])}
    byseq = defaultdict(list)
    invalid = defaultdict(int)
    for gid, row in genes.items():
        sid = _sid(row)
        if sid not in selected:
            continue
        try:
            start, end = int(row['start0']), int(row['end0'])
            length = int(seq[sid]['length'])
            assert 0 <= start < end <= length and row['strand'] in {'+', '-'}
        except (KeyError, ValueError, TypeError, AssertionError):
            invalid[sid] += 1
            continue
        byseq[sid].append(row)
    chosen = sorted([sid for sid in selected if byseq[sid] and not invalid[sid]],
                    key=lambda s: (-int(seq[s]['length']), s))[:2]
    if not chosen:
        return _empty_panel(ctx, stem, title,
                            'No selected sequence has a complete, valid gene-coordinate table for illustration.', caption)
    plt = ctx['plt']
    fig, axes = plt.subplots(len(chosen), 1, figsize=(7.08, 2.05 * len(chosen) + .65), squeeze=False)
    exported = []
    for ax, sid in zip(axes[:, 0], chosen):
        length = int(seq[sid]['length'])
        ordered = sorted(byseq[sid], key=lambda r: (int(r['start0']), int(r['end0']), r['gene_id']))
        highlights = sorted([r for r in ordered if annotations.get(r['gene_id'])],
                            key=lambda r: (0 if any(k in annotations[r['gene_id']] for k in ['Defense', 'Anti-defense', 'KO']) else 1,
                                           int(r['start0']), r['gene_id']))[:6]
        highlight_ids = {r['gene_id'] for r in highlights}
        label_order = {r['gene_id']: i for i, r in enumerate(sorted(highlights, key=lambda r: int(r['start0'])))}
        ax.axhline(0, color=GRID, linewidth=.8, zorder=0)
        for row in ordered:
            gid = row['gene_id']
            start, end = int(row['start0']), int(row['end0'])
            forward = row['strand'] == '+'
            y = .12 if forward else -.12
            parts = annotations.get(gid, {})
            color, display_type = _gene_color(parts)
            x = start if forward else end
            dx = end - start if forward else start - end
            ax.arrow(x / 1000, y, dx / 1000, 0, width=.085,
                     head_width=.20, head_length=min(length * .008, abs(dx) * .45) / 1000,
                     length_includes_head=True, facecolor=color, edgecolor='white', linewidth=.25,
                     zorder=3)
            if gid in highlight_ids:
                index = label_order[gid]
                label_y = .63 + .27 * (index % 2)
                # Spread the few labels into ordered slots rather than stacking
                # adjacent short genes' text at almost identical coordinates.
                label_x = length * (index + .5) / len(highlights) / 1000
                ax.annotate(_gene_label(parts), xy=((start + end) / 2000, y + .07),
                            xytext=(label_x, label_y), fontsize=6.5,
                            color=INK, ha='center', va='bottom',
                            arrowprops=dict(arrowstyle='-', color='#91A0AB', linewidth=.6))
            exported.append(dict(viral_sequence_id=sid, display_alias=ctx['aliases'].get(sid, sid),
                                 gene_id=gid, start0=start, end0=end, strand=row['strand'],
                                 sequence_length=length, display_type=display_type,
                                 labeled=str(gid in highlight_ids).lower(),
                                 phrog_ids=';'.join(sorted(parts.get('PHROG', []))),
                                 defense_subtypes=';'.join(sorted(parts.get('Defense', []))),
                                 anti_defense_subtypes=';'.join(sorted(parts.get('Anti-defense', []))),
                                 ko_ids=';'.join(sorted(parts.get('KO', []))),
                                 interpretation='own_sequence_computational_evidence;unreviewed_illustrative_example'))
        alias = ctx['aliases'].get(sid, sid)
        ax.set_title(f'{alias}  |  {length / 1000:.1f} kb  |  {len(ordered)} predicted genes',
                     loc='left', fontsize=8.5, color=INK, pad=8, fontweight='bold')
        ax.set_xlim(-length / 1000 * .025, length / 1000 * 1.025)
        ax.set_ylim(-.48, 1.24)
        ax.set_yticks([])
        ax.set_xlabel('Position in this viral sequence (kb)', color=MUTED, fontsize=7.5)
        ax.tick_params(axis='x', labelsize=7, colors=MUTED, length=3)
        for side in ['top', 'left', 'right']:
            ax.spines[side].set_visible(False)
        ax.spines['bottom'].set_color(GRID)
    fig.suptitle(title, x=.07, y=.99, ha='left', fontsize=11, color=INK, fontweight='bold')
    handles = [plt.Line2D([], [], marker='s', linestyle='none', markersize=5, color=c)
               for c in [COLORS[0], COLORS[1], COLORS[2], COLORS[3], '#BCC7CF']]
    fig.legend(handles, ['PHROG', 'Defense', 'Anti-defense', 'KO', 'Other called genes'],
               loc='lower center', ncol=5, frameon=False, fontsize=7,
               handletextpad=.35, columnspacing=1.15, bbox_to_anchor=(.51, .007))
    fig.subplots_adjust(left=.08, right=.96, top=.87, bottom=.14, hspace=.80)
    ctx['save'](fig, stem)
    ctx['export'](stem, exported)
    return dict(stem=stem, title=title, caption=caption,
                selection='Among the pooled-occurrence display set, choose up to two sequences by descending length, then ID; require every coordinate to be valid.')


def draw_function_figures(ctx):
    """Return display metadata; all figures use actual same-sequence evidence.

    Optional input keys ``phrog_hits``, ``systems`` and ``function_status`` are
    intentionally separate from the compact raw intersection evidence table.
    ``save(fig, stem)`` owns figure closure and SVG/PDF/PNG output.
    ``export(stem, rows)`` writes the exact machine-readable plotted evidence.
    """
    selected = list(dict.fromkeys(ctx.get('selected_m_sequences', [])))[:16]
    genes, evidence, annotations = _evidence(ctx)
    return [_dot_matrix(ctx, selected, evidence),
            _gene_tracks(ctx, selected, genes, annotations)]
