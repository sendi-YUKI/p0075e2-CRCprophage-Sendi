"""Select an output's explicit Nextflow session UUID, never CLI 'last'/run name.

Nextflow 25.10.4 may normalize ``-resume RUN_NAME`` to ``last`` (upstream #3362).
Exact history lookup is therefore wrapper-side, and postrun session equality is
mandatory. This module performs no writes and never launches Nextflow.
"""
from pathlib import Path
import re

UUID_RE = re.compile(r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}')


def checked_uuid(value):
    if not isinstance(value, str) or not UUID_RE.fullmatch(value.lower()):
        raise RuntimeError('Nextflow resume/session value must be a full UUID, never last or a run name')
    return value.lower()


def history_session(history_path, run_name):
    """Read the exact RUN_NAME column (2), not a substring in command text."""
    if not isinstance(run_name, str) or not run_name or run_name == 'last' or any(c in run_name for c in '\t\r\n'):
        raise RuntimeError('Invalid exact Nextflow run name for history lookup')
    history = Path(history_path)
    if not history.is_file():
        return None
    found = set()
    for line in history.read_text(encoding='utf-8').splitlines():
        fields = line.split('\t', 6)
        if len(fields) >= 3 and fields[2] == run_name:
            if len(fields) < 7:
                raise RuntimeError('Malformed matching Nextflow history record')
            found.add(checked_uuid(fields[5]))
    if len(found) > 1:
        raise RuntimeError('Ambiguous Nextflow history: exact run name has multiple session UUIDs')
    return next(iter(found)) if found else None


def select_resume_session(previous, workflow_status, history_path):
    """Resolve before replacing state; keep the UUID across prelaunch failures."""
    pinned = previous.get('resume_session_id')
    if pinned is not None:
        pinned = checked_uuid(pinned)
    observed = previous.get('nextflow_session_id')
    if observed is not None:
        observed = checked_uuid(observed)
    name = previous.get('nextflow_run_name') if previous.get('nextflow_started') else None
    if name:
        receipt = None
        if workflow_status.get('run_name') == name:
            receipt = checked_uuid(workflow_status.get('session_id'))
        from_history = history_session(history_path, name)
        evidence = {x for x in (observed, receipt, from_history) if x}
        if len(evidence) > 1:
            raise RuntimeError('Conflicting Nextflow session evidence for this output run')
        if evidence:
            result = evidence.pop()
            if pinned and result != pinned:
                raise RuntimeError('Recorded Nextflow run used a different session than its requested resume UUID')
            return result
        # A wrapper can fail just before Nextflow records its new run name.
        # The previously saved explicit UUID remains the only permissible target.
        if pinned:
            return pinned
        raise RuntimeError('Cannot identify this output run in its receipt or exact Nextflow history')
    if pinned or observed:
        if pinned and observed and pinned != observed:
            raise RuntimeError('Conflicting saved Nextflow UUIDs')
        return pinned or observed
    legacy = previous.get('resume_target')
    if legacy:
        if isinstance(legacy, str) and UUID_RE.fullmatch(legacy.lower()):
            return checked_uuid(legacy)
        target = history_session(history_path, legacy)
        if workflow_status.get('run_name') == legacy:
            receipt = checked_uuid(workflow_status.get('session_id'))
            if target and target != receipt:
                raise RuntimeError('Conflicting legacy Nextflow history and receipt')
            return receipt
        if target:
            return target
        raise RuntimeError('Cannot resolve saved legacy run name to an exact session UUID')
    # No prior scientific launch: a failed first preflight may safely start its
    # first session, but an unbound workflow receipt must not be borrowed.
    if workflow_status:
        raise RuntimeError('Workflow receipt has no matching previous output run identity')
    return None


def verify_nextflow_session(workflow_status, run_name, requested_uuid=None):
    if workflow_status.get('run_name') != run_name:
        raise RuntimeError('Workflow receipt does not belong to this Nextflow run name')
    actual = checked_uuid(workflow_status.get('session_id'))
    if requested_uuid is not None and actual != checked_uuid(requested_uuid):
        raise RuntimeError('Nextflow resumed a different session UUID than requested')
    return actual
