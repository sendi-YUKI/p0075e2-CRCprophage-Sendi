#!/usr/bin/env python3
"""Turn the protocol's banned legacy behaviours into assertions that fail.

Protocol section 12 lists six things the new pipeline must not do. Until now
they existed only as prose, while the behaviours themselves are hard-coded
first-class grouping variables in old_version/r_code. Section 12's stated plan
is to reuse the existing pipeline rather than rebuild it, so reuse is a given
and prose is not a control.

Each check below fails loudly rather than warning. Run it in CI and before any
freeze.

Usage:
    python3 bin/legacy_guard.py [--repo PATH] [--json OUT]
Exit code 0 when clean, 1 when any ban is violated.
"""
import argparse
import json
import re
import sys
from pathlib import Path

# Files that are historical evidence, not executable pipeline code.
EXCLUDE_DIRS = {'.git', '__pycache__', 'envs', 'results', 'old_version',
                'baseline_evidence', 'data', '.pixi'}
CODE_SUFFIX = {'.py', '.R', '.r', '.nf', '.config', '.sh'}


def code_files(repo):
    for p in repo.rglob('*'):
        if p.suffix not in CODE_SUFFIX or not p.is_file():
            continue
        if any(part in EXCLUDE_DIRS for part in p.parts):
            continue
        yield p


def scan(repo, pattern, flags=re.I):
    rx = re.compile(pattern, flags)
    hits = []
    for p in code_files(repo):
        try:
            text = p.read_text(errors='replace')
        except OSError:
            continue
        for n, line in enumerate(text.splitlines(), 1):
            if line.lstrip().startswith(('#', '//', '--')):
                continue          # a comment naming the ban is not a violation
            if rx.search(line):
                hits.append('%s:%d: %s' % (p.relative_to(repo), n, line.strip()[:120]))
    return hits


# Explicit exemptions. Each entry must name the check, the exact path and why.
# Nothing is exempt by pattern or by looking like a test; an exempt hit is still
# printed, it just does not fail the run.
EXEMPT = [
    dict(check='ban1_species_to_role', path='bin/legacy_guard.py',
         why='the guard itself has to name the banned words to search for them'),
    dict(check='ban4_breadth10_single_mapq', path='bin/legacy_guard.py',
         why='the guard itself has to name the banned threshold to search for it'),
    dict(check='ban4_breadth10_single_mapq', path='bin/run_stage.py',
         why='docstring of research_profile_guard explains the ban it enforces'),
    dict(check='ban4_breadth10_single_mapq',
         path='bootstrap/spark/basic_build_community_smoke.py',
         why=('synthetic engineering smoke fixture, status=synthetic, not a research path; '
              'research runs are gated by --research-run and research_eligible')),
]


def exempt_reason(check_id, relpath):
    for e in EXEMPT:
        if e['check'] == check_id and e['path'] == relpath:
            return e['why']
    return None


CHECKS = [
    dict(
        id='ban1_species_to_role',
        ban='不能把某 species 自动标为 oral / protector',
        why=('物种到 role 的三分类是旧分析的一级分组变量'
             '（old_version/r_code/aim_beta.R:978,988），oral/gut 是正文图的分面变量'
             '（aim_beta_figures.R:173）。role 只能来自 metadata_evidence 且证据级别 A 或 B。'),
        pattern=r'\b(oncomicrobe|protective|protector)\b',
    ),
    dict(
        id='ban2_legacy_626_branch',
        ban='不能用旧 626 分支的功能补新目录',
        why='所有主分析实例与功能必须来自同一个 master universe。',
        pattern=r'\b626[_-]?(genome|branch|amg|cargo)\b',
    ),
    dict(
        id='ban3_species_ani_dereplication',
        ban='不能按 species 做 95% ANI 去重而删掉种内多样性',
        why='种内变异正是本项目的分析单位；去重会删掉研究对象本身。',
        pattern=r'\bderepl?icat\w*|\bdrep\b',
    ),
    dict(
        id='ban4_breadth10_single_mapq',
        ban='不能用 10% reference breadth 加单个 MAPQ 阈值宣称可靠的整条 phage 检出',
        why=('检出判定必须同时读 breadth、identity、MAPQ 与可区分区域长度四个键。'
             'dev profile 保留此行为但必须带 research_eligible: false。'),
        pattern=r'breadth\w*\s*[=:<>]\s*0?\.1\b(?!\d)',
    ),
    dict(
        id='ban5_nm_over_seqlen_identity',
        ban='不能把局部比对里的 1 − NM/len(SEQ) 当作无需校验的 identity',
        why='local 模式下用完整 SEQ 长度做分母会稀释错配；必须基于有效 CIGAR 与 NM 并说明 indel 与 clip 如何计入。',
        pattern=r'1\s*-\s*NM\s*/\s*(len\s*\(\s*SEQ|seq_len|length\s*\(\s*seq)',
    ),
    dict(
        id='ban6_contig_remainder_as_host_backbone',
        ban='不能用包含其他 prophage 的同 contig 余区直接称纯宿主 backbone',
        why='宿主 backbone 必须从同 contig 扣除全部已识别移动区，并输出扣除清单。',
        pattern=r'(backbone|host_core)\w*\s*=.*(contig|remainder|flank)(?!.*mask)',
    ),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument('--json')
    a = ap.parse_args()
    repo = Path(a.repo).resolve()

    report, failed = [], 0
    for c in CHECKS:
        raw = scan(repo, c['pattern'])
        hits, exempted = [], []
        for h in raw:
            rel = h.split(':', 1)[0]
            why = exempt_reason(c['id'], rel)
            (exempted if why else hits).append(h if not why else '%s  [豁免: %s]' % (h, why))
        rec = dict(id=c['id'], ban=c['ban'], why=c['why'], pattern=c['pattern'],
                   n_hits=len(hits), hits=hits[:20],
                   n_exempt=len(exempted), exempt=exempted[:20])
        report.append(rec)
        for e in exempted:
            print('     豁免 %s' % e)
        if hits:
            failed += 1
            print('FAIL %s' % c['id'])
            print('  禁令: %s' % c['ban'])
            print('  理由: %s' % c['why'])
            for h in hits[:10]:
                print('    %s' % h)
            if len(hits) > 10:
                print('    ... 另有 %d 处' % (len(hits) - 10))
        else:
            print('ok   %s' % c['id'])

    if a.json:
        Path(a.json).write_text(json.dumps(
            dict(repo=str(repo), n_checks=len(CHECKS), n_failed=failed, checks=report),
            indent=2, ensure_ascii=False) + '\n')
    print('\n%d/%d 条禁令被违反' % (failed, len(CHECKS)))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
