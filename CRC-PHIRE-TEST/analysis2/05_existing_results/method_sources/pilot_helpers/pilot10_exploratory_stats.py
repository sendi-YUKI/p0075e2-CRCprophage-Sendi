#!/usr/bin/env python3
"""Fixed ten-donor exploratory tests; standard library, no upstream recomputation.

Production requires Slurm. P/Q values describe this discovery cohort only: the
label enumeration is a small-sample reference calculation, not a randomized
experiment or a way to remove confounding. Each donor is counted once per vOTU.
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import html
import itertools
import json
import math
import os
from pathlib import Path
import re
import statistics
import csv

DEFAULT_ROOT = Path('/srv/CRC-PHIRE/analysis1')
EVIDENCE_FIELDS = ['has_reference_predicted_provirus', 'has_metagenome_predicted_provirus',
                  'has_metagenome_two_sided_host_markers', 'evidence_label',
                  'n_reference_provirus_occurrences', 'n_metagenome_provirus_occurrences']
DIFFERENTIAL_FIELDS = [
    'votu_id', 'reference_id', 'crc_detection_n', 'control_detection_n',
    'crc_detected', 'control_detected', 'crc_detection_fraction', 'control_detection_fraction',
    'detection_fraction_difference', 'fisher_p', 'fisher_q_global', 'fisher_q_endpoint',
    'detection_test_status', 'crc_abundance_n', 'control_abundance_n',
    'crc_median_rpkm', 'control_median_rpkm', 'median_rpkm_difference', 'cliffs_delta',
    'abundance_p', 'abundance_q_global', 'abundance_q_endpoint', 'abundance_test_status',
    'abundance_permutations', *EVIDENCE_FIELDS]
DONOR_FIELDS = ['votu_id', 'reference_id', 'biological_sample_id', 'group',
                'detection_value', 'abundance_rpkm', 'detection_valid', 'abundance_valid',
                'detection_missing_reason', 'abundance_missing_reason',
                'measurement_status', 'detection_status']
FUNCTION_FIELDS = [
    'votu_id', 'reference_id', 'viral_sequence_id', 'membership_status',
    'representative_id', 'annotation_type', 'feature', 'gene_ids',
    'evidence_file', 'evidence_sha256', 'sequence_source_sides',
    'metagenome_discovery_donors', 'fisher_p', 'fisher_q_global',
    'abundance_p', 'abundance_q_global', 'cliffs_delta', 'evidence_label', 'interpretation']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def table(path, fields=()):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(set(fields) <= set(reader.fieldnames or []), f'Missing columns: {path}: {fields}')
        return list(reader)


def write_table(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def verify_outputs(base, mapping, inputs):
    require(isinstance(mapping, dict) and bool(mapping), f'Missing output hashes under {base}')
    base = Path(base).resolve()
    for rel, digest in mapping.items():
        path = (base / rel).resolve()
        require(path.is_relative_to(base), f'Unsafe receipt path: {rel}')
        require(path.is_file() and sha(path) == digest, f'Hash mismatch: {path}')
        inputs[str(path)] = digest


def fisher_exact_two_sided(a, b, c, d):
    """Probability-ordering two-sided Fisher test using exact integer weights."""
    cells = (a, b, c, d)
    require(all(isinstance(x, int) and x >= 0 for x in cells), 'Invalid 2x2 table')
    n1, n2, total_success = a + b, c + d, a + c
    require(n1 > 0 and n2 > 0, 'Both groups need observations')
    denominator = math.comb(n1 + n2, total_success)
    observed_weight = math.comb(n1, a) * math.comb(n2, c)
    numerator = 0
    for x in range(max(0, total_success - n2), min(n1, total_success) + 1):
        weight = math.comb(n1, x) * math.comb(n2, total_success - x)
        if weight <= observed_weight:
            numerator += weight
    return min(1.0, numerator / denominator)


def doubled_average_ranks(values):
    """Integer midranks avoid floating equality ambiguity when ties occur."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0] * len(values)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and values[order[stop]] == values[order[start]]:
            stop += 1
        for i in order[start:stop]:
            ranks[i] = start + stop + 1
        start = stop
    return ranks


def exact_rank_permutation(x, y):
    """Enumerate all group assignments, including the observed one; ties valid.

    The absolute difference of average ranks is equivalent to centered U for
    fixed group sizes. No Monte Carlo sampling, pseudocount, or asymptotic P.
    """
    require(bool(x) and bool(y), 'Empty rank-test group')
    values = list(x) + list(y)
    require(all(math.isfinite(v) and v >= 0 for v in values), 'Invalid abundance')
    require(len(values) <= 10, 'This pilot permits at most ten donors')
    ranks = doubled_average_ranks(values)
    n1, n = len(x), len(values)
    total = sum(ranks)
    observed = abs(n * sum(ranks[:n1]) - n1 * total)
    extreme = 0
    count = 0
    for indices in itertools.combinations(range(n), n1):
        statistic = abs(n * sum(ranks[i] for i in indices) - n1 * total)
        extreme += statistic >= observed
        count += 1
    return extreme / count, count


def cliffs_delta(x, y):
    require(bool(x) and bool(y), 'Empty effect-size group')
    return sum((a > b) - (a < b) for a in x for b in y) / (len(x) * len(y))


def bh_adjust(p_values):
    """BH over supplied valid endpoints; caller keeps unavailable tests as NA."""
    require(all(math.isfinite(p) and 0 <= p <= 1 for p in p_values), 'Invalid P value')
    ordered = sorted(range(len(p_values)), key=lambda i: p_values[i])
    result = [None] * len(p_values)
    running = 1.0
    for rank0 in range(len(ordered) - 1, -1, -1):
        index = ordered[rank0]
        running = min(running, p_values[index] * len(ordered) / (rank0 + 1))
        result[index] = min(1.0, running)
    return result


def finite_nonnegative(value):
    if value is None or str(value).strip().lower() in {'', 'na', 'nan', 'none', 'null'}:
        return None, 'missing_abundance'
    try:
        result = float(value)
    except (ValueError, TypeError):
        return None, 'invalid_abundance_number'
    if not math.isfinite(result) or result < 0:
        return None, 'abundance_not_finite_nonnegative'
    return result, ''


def donor_values(row):
    assessed = row['measurement_status'] == 'assessed'
    detection = row['detection_status']
    require(assessed or detection not in {'detected', 'not_detected'},
            'Not-assessed sample has an evaluated detection')
    valid_detection = assessed and detection in {'detected', 'not_detected'}
    abundance, reason = finite_nonnegative(row.get('fragment_rpkm_all_clean'))
    if not assessed:
        abundance, reason = None, 'measurement_not_assessed'
    return dict(votu_id=row['votu_id'], reference_id=row['reference_id'],
                biological_sample_id=row['biological_sample_id'], group=row['group'],
                detection_value=int(detection == 'detected') if valid_detection else None,
                abundance_rpkm=abundance,
                detection_valid=str(valid_detection).lower(),
                abundance_valid=str(abundance is not None).lower(),
                detection_missing_reason='' if valid_detection else
                    ('measurement_not_assessed' if not assessed else 'detection_status_' + detection),
                abundance_missing_reason=reason, measurement_status=row['measurement_status'],
                detection_status=detection)


def analyze(rows, groups, evidence):
    require(len(groups) == 10 and Counter(groups.values()) == {'CRC': 5, 'control': 5},
            'Expected frozen ten donors, five per group')
    byvotu, targets, unique = defaultdict(list), {}, set()
    for source in rows:
        vid, target, donor = source['votu_id'], source['reference_id'], source['biological_sample_id']
        require(bool(vid) and bool(target), 'Blank measured vOTU/target')
        require(donor in groups and source['group'] == groups[donor], 'Frozen group mismatch')
        require((vid, donor) not in unique, 'Duplicate donor-vOTU (pseudoreplication)')
        unique.add((vid, donor))
        require(vid not in targets or targets[vid] == target, 'One vOTU maps to multiple targets')
        targets[vid] = target
        byvotu[vid].append(donor_values(source))
    require(len(set(targets.values())) == len(targets), 'One target maps to multiple vOTUs')
    require(set(byvotu) <= set(evidence), 'Missing evidence rows for measured vOTUs')
    output, donors = [], []
    for vid, records in sorted(byvotu.items()):
        require(len(records) == 10 and {r['biological_sample_id'] for r in records} == set(groups),
                f'Incomplete ten-donor rectangle: {vid}; absent input row is not a negative')
        current = dict.fromkeys(DIFFERENTIAL_FIELDS)
        current.update(votu_id=vid, reference_id=targets[vid])
        current.update({field: evidence[vid][field] for field in EVIDENCE_FIELDS})
        for prefix, group in [('crc', 'CRC'), ('control', 'control')]:
            subset = [r for r in records if r['group'] == group]
            detections = [r['detection_value'] for r in subset if r['detection_valid'] == 'true']
            abundances = [r['abundance_rpkm'] for r in subset if r['abundance_valid'] == 'true']
            current[prefix + '_detection_n'] = len(detections)
            current[prefix + '_detected'] = sum(detections)
            current[prefix + '_detection_fraction'] = statistics.mean(detections) if detections else None
            current[prefix + '_abundance_n'] = len(abundances)
            current[prefix + '_median_rpkm'] = statistics.median(abundances) if abundances else None
        a, b = current['crc_detection_fraction'], current['control_detection_fraction']
        current['detection_fraction_difference'] = a - b if a is not None and b is not None else None
        if min(current['crc_detection_n'], current['control_detection_n']) >= 3:
            a, c = current['crc_detected'], current['control_detected']
            current['fisher_p'] = fisher_exact_two_sided(
                a, current['crc_detection_n'] - a, c, current['control_detection_n'] - c)
            current['detection_test_status'] = 'tested_fisher_two_sided'
        else:
            current['detection_test_status'] = 'not_tested_fewer_than_3_evaluable_per_group'
        x = [r['abundance_rpkm'] for r in records if r['group'] == 'CRC' and r['abundance_valid'] == 'true']
        y = [r['abundance_rpkm'] for r in records if r['group'] == 'control' and r['abundance_valid'] == 'true']
        if x and y:
            current['median_rpkm_difference'] = statistics.median(x) - statistics.median(y)
            current['cliffs_delta'] = cliffs_delta(x, y)
        if min(len(x), len(y)) >= 3:
            current['abundance_p'], current['abundance_permutations'] = exact_rank_permutation(x, y)
            current['abundance_test_status'] = ('tested_all_values_equal_p1' if len(set(x + y)) == 1
                                                else 'tested_exact_two_sided_rank_permutation')
        else:
            current['abundance_test_status'] = 'not_tested_fewer_than_3_valid_per_group'
        output.append(current)
        donors.extend(sorted(records, key=lambda r: r['biological_sample_id']))
    # One primary family includes both endpoints across every testable vOTU.
    tests = [(r, endpoint) for r in output for endpoint in ('fisher', 'abundance')
             if r[endpoint + '_p'] is not None]
    for (row, endpoint), q in zip(tests, bh_adjust([r[e + '_p'] for r, e in tests])):
        row[endpoint + '_q_global'] = q
    for endpoint in ('fisher', 'abundance'):
        valid = [r for r in output if r[endpoint + '_p'] is not None]
        for row, q in zip(valid, bh_adjust([r[endpoint + '_p'] for r in valid])):
            row[endpoint + '_q_endpoint'] = q
    return output, donors


def function_links(catalog_root, differential, inputs):
    """Join an own-sequence annotation to its assigned member; never broadcast."""
    paths = {name: catalog_root / rel for name, rel in {
        'members': 'catalog/vOTU_members.tsv', 'occurrences': 'catalog/source_occurrences.tsv',
        'genes': 'genes/gene_table.tsv', 'phrog': 'functions/function_hits.tsv',
        'systems': 'functions/systems.tsv',
        'ko': 'viral_kofam/kofam/viral_member_ko_counts.tsv'}.items()}
    for path in paths.values():
        require(path.is_file() and str(path.resolve()) in inputs,
                f'Functional input not covered by terminal receipt: {path}')
    differential = {r['votu_id']: r for r in differential}
    members = {}
    for row in table(paths['members'], ['viral_sequence_id', 'votu_id']):
        sid = row['viral_sequence_id']
        require(sid not in members, f'Duplicate member sequence: {sid}')
        members[sid] = row
    genes = {}
    for row in table(paths['genes'], ['gene_id', 'viral_sequence_id']):
        require(row['gene_id'] not in genes, 'Duplicate gene ID')
        require(row['viral_sequence_id'] in members, 'Unknown gene source sequence')
        genes[row['gene_id']] = row['viral_sequence_id']
    sources, discovery = defaultdict(set), defaultdict(set)
    for row in table(paths['occurrences'], ['viral_sequence_id', 'source_kind', 'biological_sample_id']):
        sid = row['viral_sequence_id']
        side = 'reference' if row['source_kind'] == 'legacy_bacterial_genome' else 'metagenome'
        sources[sid].add(side)
        if side == 'metagenome':
            discovery[sid].add(row['biological_sample_id'])
    result = []

    def add(sid, annotation_type, feature, gene_ids, source, meaning):
        require(sid in members, f'Unknown annotated sequence: {sid}')
        member = members[sid]
        vid = member['votu_id']
        for gid in re.split(r'[;,|\s]+', gene_ids or ''):
            if gid:
                require(gid in genes and genes[gid] == sid, 'Annotation gene/sequence mismatch')
        if vid not in differential:
            return
        stats = differential[vid]
        result.append(dict(votu_id=vid, reference_id=stats['reference_id'], viral_sequence_id=sid,
                           membership_status=member.get('membership_status', ''),
                           representative_id=member.get('representative_id', ''),
                           annotation_type=annotation_type, feature=feature, gene_ids=gene_ids,
                           evidence_file=str(source), evidence_sha256=inputs[str(source.resolve())],
                           sequence_source_sides=';'.join(sorted(sources[sid])),
                           metagenome_discovery_donors=';'.join(sorted(discovery[sid])),
                           **{key: stats[key] for key in ['fisher_p', 'fisher_q_global', 'abundance_p',
                                                       'abundance_q_global', 'cliffs_delta', 'evidence_label']},
                           interpretation=meaning + ';vOTU_test_not_gene_carriage_or_expression_test;'
                               'discovery_donors_not_all_gene_carriers;no_annotation_transfer'))
    for row in table(paths['phrog']):
        if row.get('accepted') == 'true' and row.get('is_best') == 'true':
            require(row['gene_id'] in genes, 'Unknown PHROG gene')
            add(genes[row['gene_id']], 'PHROG_homology', row.get('phrog', ''), row['gene_id'],
                paths['phrog'], 'own_sequence_accepted_best_homology_not_experimental_function')
    for row in table(paths['systems']):
        sid = row.get('viral_sequence_id') or row.get('contig_id')
        add(sid, row.get('activity', 'system_prediction'), row.get('subtype') or row.get('type', ''),
            row.get('all_gene_ids') or row.get('required_gene_ids') or '', paths['systems'],
            'own_sequence_system_prediction_not_validated_defense_or_antidefense_activity')
    for row in table(paths['ko']):
        sid = row['viral_sequence_id']
        require(sid in members and row.get('votu_id', '') == members[sid]['votu_id'], 'KO member mismatch')
        add(sid, 'KO_homology', row['ko'], row.get('gene_ids', ''), paths['ko'],
            'own_sequence_KO_homology_candidate_cargo_not_confirmed_AMG')
    return sorted(result, key=lambda r: (r['votu_id'], r['viral_sequence_id'], r['annotation_type'], r['feature'], r['gene_ids']))


def run(root, *, _synthetic_test=False):
    require(bool(os.environ.get('SLURM_JOB_ID')), 'Production execution requires SLURM_JOB_ID')
    root = Path(root).resolve()
    report = root / '12_reports/pilot10'
    catroot = root / '06_viral_catalog/pilot10_combined'
    evidence_dir = report / 'prophage_evidence'
    final_path, evidence_path = report / 'pilot10_final_acceptance.json', evidence_dir / 'evidence_manifest.json'
    frozen_path = root / '01_manifests/pilot10_frozen_assemblies.json'
    inputs = {str(path): sha(path) for path in [final_path, evidence_path, frozen_path, Path(__file__).resolve()]}
    final, evidence_receipt, frozen = read(final_path), read(evidence_path), read(frozen_path)
    require(bool(final.get('synthetic_preview', False)) == _synthetic_test,
            'Synthetic fixture and production data must not be mixed')
    if _synthetic_test:
        require((root / 'SYNTHETIC_TEST_DATA_ONLY.txt').is_file(), 'Missing synthetic fixture marker')
    require(final['status'] == 'completed_balanced_engineering_pilot' and final['donors'] == 10,
            'Ten-donor report not accepted')
    verify_outputs(report, final['output_sha256'], inputs)
    require('donor_votu_measurement.tsv' in final['output_sha256'], 'Donor measurements not receipted')
    require(evidence_receipt['status'] == 'completed', 'Prophage evidence not accepted')
    verify_outputs(evidence_dir, evidence_receipt['outputs'], inputs)
    require('votu_evidence.tsv' in evidence_receipt['outputs'], 'vOTU evidence not receipted')
    for source, digest in evidence_receipt['input_sha256'].items():
        path = Path(source)
        if not path.is_absolute():
            path = root / path
        require(path.is_file() and sha(path) == digest, f'Evidence source changed: {path}')
        inputs[str(path.resolve())] = digest
    stage, terminal = catroot / 'pipeline_info/stage_state.json', catroot / 'pipeline_info/terminal_output_manifest.json'
    from pilot10_catalog_acceptance import verified_catalog_completion
    publication = verified_catalog_completion(root, catroot)
    if publication is not None:
        inputs[str(publication)] = sha(publication)
    inputs.update({str(stage): sha(stage), str(terminal): sha(terminal)})
    verify_outputs(catroot, read(terminal), inputs)
    for relative, hash_field in [('genes/gene_calling_status.json', 'output_files'),
                                 ('functions/viral_function_status.json', 'output_files')]:
        receipt_path = catroot / relative
        require(str(receipt_path) in inputs, f'Child receipt not in terminal manifest: {relative}')
        receipt = read(receipt_path)
        require(receipt['status'] == 'completed', f'Incomplete child artifact: {relative}')
        verify_outputs(receipt_path.parent, receipt[hash_field], inputs)
    # run_stage currently adds this optional artifact after terminal validation.
    # Also support receipts made before that addition, anchored to the exact
    # already-verified genes and membership rather than a filename assumption.
    ko_receipt_path = catroot / 'viral_kofam/kofam/viral_ko_manifest.json'
    ko_receipt = read(ko_receipt_path)
    require(ko_receipt['status'] == 'completed', 'Incomplete KOfam artifact')
    for field, relative in [('genes_sha256', 'genes/gene_table.tsv'),
                            ('membership_sha256', 'catalog/vOTU_members.tsv')]:
        anchored_path = catroot / relative
        require(str(anchored_path) in inputs, f'Missing KOfam input anchor: {relative}')
        require(ko_receipt[field] == inputs[str(anchored_path)], f'KOfam used different {relative}')
    inputs[str(ko_receipt_path)] = sha(ko_receipt_path)
    verify_outputs(ko_receipt_path.parent, ko_receipt['outputs'], inputs)
    evidence = {}
    for row in table(evidence_dir / 'votu_evidence.tsv', ['votu_id', *EVIDENCE_FIELDS]):
        require(row['votu_id'] not in evidence, 'Duplicate evidence vOTU')
        for field in EVIDENCE_FIELDS[:3]:
            require(row[field] in {'true', 'false', 'unknown'}, f'Invalid evidence state: {field}')
        evidence[row['votu_id']] = row
    raw = table(report / 'donor_votu_measurement.tsv',
                ['votu_id', 'reference_id', 'biological_sample_id', 'group', 'measurement_status', 'detection_status'])
    differential, donors = analyze(raw, frozen['groups'], evidence)
    function_rows = function_links(catroot, differential, inputs)
    destination = report / 'exploratory'
    if destination.exists():
        old_path = destination / 'exploratory_manifest.json'
        require(old_path.is_file(), 'Partial exploratory output exists; preserve and investigate')
        old = read(old_path)
        require(old['status'] == 'completed' and old['input_sha256'] == inputs,
                'Existing exploratory result has different inputs/code; preserve it')
        verify_outputs(destination, old['outputs'], {})
        print(json.dumps({'status': 'verified_existing_exploratory_result', 'path': str(destination)}))
        return old
    staging = report / ('.exploratory_build_' + str(os.environ['SLURM_JOB_ID']))
    require(not staging.exists(), f'Existing partial staging output; preserve: {staging}')
    staging.mkdir()
    write_table(staging / 'votu_differential.tsv', differential, DIFFERENTIAL_FIELDS)
    write_table(staging / 'donor_analysis_values.tsv', donors, DONOR_FIELDS)
    write_table(staging / 'differential_function_links.tsv', function_rows, FUNCTION_FIELDS)
    tests = sum(row[e + '_p'] is not None for row in differential for e in ('fisher', 'abundance'))
    discoveries = sum(row[e + '_q_global'] is not None and row[e + '_q_global'] <= .05
                      for row in differential for e in ('fisher', 'abundance'))
    md = '\n\n'.join([
        '# CRC-PHIRE：10人探索性差异分析',
        'SYNTHETIC SOFTWARE TEST ONLY — NOT BIOLOGICAL RESULTS.' if _synthetic_test else '',
        '本分析补充比较Feng队列5例CRC与5名研究对照；十人同时参与病毒发现与测量，属于小样本观察性先导。'
        '对照并非全部健康人。结果不能代替宿主校正、独立队列验证或机制实验。',
        f'保留全部{len(differential)}个已测量vOTU，共{len(donors)}条供者×vOTU记录；'
        f'共有{tests}个可计算端点检验，主全局BH Q≤0.05的端点数为{discoveries}。零个通过同样是完整有效的输出，不能据此证明两组相同。',
        '## 统计单位与缺失值',
        '每位供者在每个vOTU中只出现一次，PE、SE、孤儿reads及多个组装实例均不作为独立患者。'
        '原始十人及分组冻结不变；未评价、缺失、无效丰度保留NA及原因，不填零。每个端点两组各至少3个有效值才检验；'
        '不满足仍保留效应摘要和不检验原因。真零值参与检验，全部同值时P=1。未检出不等于生物学不存在。',
        '## 两个预设端点',
        '检出比例使用双侧Fisher精确检验，报告检出人数、有效分母及CRC−对照比例差。'
        '丰度使用fragment_rpkm_all_clean，以平均秩差的绝对值为统计量，完全枚举所有标签分配（最多252种）；'
        '允许并列秩，等价于固定组大小下的中心化Mann–Whitney U，无渐近近似或随机抽样。'
        '报告两组中位数、中位数差及Cliff’s delta；正值表示本先导CRC侧更高，无任意伪计数或无限倍数变化。'
        '丰度用所有可评价的连续RPKM（包括未达到广泛检出阈值的非零值），不会因为检出标记为阴性强制置零。',
        '标签置换的精确枚举只表示该统计量在交换性零假设下的计算精确；患者不是随机分组，'
        '年龄、共病、药物、谱系、宿主丰度和其他混杂仍可产生差异。P值不证明CRC因果作用。',
        '## 多重检验与prophage分层',
        '主结果为fisher_q_global与abundance_q_global：全部vOTU、两个可计算端点的P值共同做一次BH校正。'
        '同时保留端点内BH列作为次要审阅字段，不能选择更小Q值替代主结果。BH的FDR解释仍依赖相应检验和依赖结构条件；'
        '病毒特征相关、样本少与目录由同批样本发现使这些Q值只能用于先导候选排序。'
        '按参考provirus、患者provirus及双侧宿主marker等证据标签分层展示，标签不改变检验集合，也不对子集另算更有利Q值。',
        '标签来自目录内至少一个位点，不能广播为每个检出该vOTU患者的整合状态。双侧marker仍是计算证据，'
        '不是独立验证；完整性标签不证明可诱导或感染能力。原始位点、边界、序列及质量证据请查prophage_evidence。',
        '## 功能怎样连接',
        'differential_function_links.tsv保留所有可链接的自身序列PHROG最佳合格命中、防御／反防御系统预测和KO同源证据，'
        '并附vOTU的统计结果；不只输出P小的对象。连接经过具体基因→病毒序列→vOTU成员关系，未将代表序列注释转移到成员。'
        '相同vOTU内基因可能不同；vOTU丰度差异不等于功能基因携带或表达差异，组装发现该序列的供者也不等于全部携带该基因的供者。'
        'KO只构成候选cargo证据，不能直接称为已确认AMG；需要结合边界和细菌侧翼排查污染。',
        '## 如何使用这次结果',
        '先看每人数据、有效样本量、效应大小和证据类型，再看主Q值。保留未显著对象、NA原因及完整表。'
        '有较大差异可以支持后续扩大队列的候选选择；无显著差异常见于样本少和多重检验，也可能是真正效应弱。'
        '本模块不拟合多协变量模型、不声称宿主校正、不生成疾病分类性能，不以先导阈值冻结正式研究结论。',
        '## 文件',
        '- votu_differential.tsv：全部vOTU的效应、P、Q、有效分母与证据分层。\n'
        '- donor_analysis_values.tsv：逐供者真实输入值、零值、缺失与原因。\n'
        '- differential_function_links.tsv：自身序列功能证据与vOTU结果连接。\n'
        '- exploratory_manifest.json：上游与输出SHA256、脚本、参数和Slurm编号。',
    ]) + '\n'
    (staging / '探索性差异分析说明.md').write_text(md, encoding='utf-8')
    page = ('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>CRC-PHIRE 10人探索性比较</title><style>body{max-width:1040px;margin:36px auto;padding:0 24px;'
            'font:17px/1.9 system-ui;color:#193645}article{white-space:pre-wrap;overflow-wrap:anywhere}</style>'
            '<h1>CRC-PHIRE 10人探索性比较</h1><article>' + html.escape(md) + '</article></html>')
    (staging / '探索性差异分析说明.html').write_text(page, encoding='utf-8')
    manifest = dict(status='completed', analysis_status='completed_exploratory_5_vs_5', created_utc=datetime.now(timezone.utc).isoformat(),
                    synthetic_preview=_synthetic_test,
                    slurm_job_id=os.environ['SLURM_JOB_ID'], code_sha256=sha(Path(__file__)),
                    donors=10, groups=dict(Counter(frozen['groups'].values())), votus=len(differential),
                    donor_rows=len(donors), function_links=len(function_rows), global_test_family_size=tests,
                    endpoints_global_q_le_0_05=discoveries, input_sha256=inputs,
                    parameters=dict(min_valid_per_group=3, detection_test='fisher_probability_ordering_two_sided',
                        abundance_test='exact_all_label_assignments_absolute_average_rank_difference',
                        abundance='fragment_rpkm_all_clean', max_permutations=252, pseudocount=None,
                        multiplicity_primary='BH_all_testable_vOTUs_both_endpoints_together',
                        multiplicity_secondary='BH_within_endpoint_not_primary',
                        evidence_subset_retesting=False, missing_imputed=False, host_adjustment=False,
                        independent_validation=False, confidence_intervals='not_estimated',
                        filter='all_measured_vOTUs_retained_no_disease_direction_filter'),
                    outputs={p.name: sha(p) for p in staging.iterdir() if p.is_file()})
    (staging / 'exploratory_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding='utf-8')
    staging.replace(destination)
    print(json.dumps({'status': manifest['status'], 'votus': len(differential), 'tests': tests,
                      'output': str(destination)}, ensure_ascii=False))
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()
    run(args.root)
