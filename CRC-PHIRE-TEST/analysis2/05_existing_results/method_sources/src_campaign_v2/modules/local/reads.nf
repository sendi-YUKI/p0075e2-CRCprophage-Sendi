process READS_CLEAN {
    tag "$id"
    input:
    tuple val(id), path(unit), path(reads)
    path profile
    path db_manifest
    path software_manifest
    path code
    val code_digest
    output:
    path "${id}", emit: units
    script:
    def human = params.human_index ? "--human-index '${params.human_index}'" : ''
    def virome = params.virome_db ? "--virome-db '${params.virome_db}'" : ''
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    # code_sha256=${code_digest}
    python3 ${code}/reads_process.py unit --unit '${unit}' --profile '${profile}' --outdir '${id}' --threads ${task.cpus} ${human} ${virome} --reference-manifest '${db_manifest}' --software-manifest '${software_manifest}'
    """
}
process READS_COLLECT {
    publishDir params.outdir, mode: 'copy'
    input:
    path units
    path profile
    path raw_manifest
    path code
    val code_digest
    output:
    path 'clean', emit: clean
    script:
    def dirs = (units instanceof List ? units : [units]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    # code_sha256=${code_digest}
    python3 ${code}/reads_process.py collect --dirs ${dirs} --outdir clean --profile '${profile}' --raw-manifest '${raw_manifest}'
    # code_sha256=${code_digest}
    python3 ${code}/reads_publish.py clean '${file(params.outdir).toAbsolutePath().normalize()}/clean'
    """
}
process READS_ASSEMBLE {
    tag "$id"
    input:
    tuple val(id), val(material), val(unit_b64), path(reads, stageAs: 'reads??/*')
    path code
    val code_digest
    output:
    tuple val(id), val(material), path("${id}"), emit: units
    script:
    def stagedReads = (reads instanceof List ? reads : [reads]).collect { it.toString() }
    def readsB64 = groovy.json.JsonOutput.toJson(stagedReads).bytes.encodeBase64().toString()
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    # code_sha256=${code_digest}
    python3 - <<'PY'
import base64, json
unit = json.loads(base64.b64decode('${unit_b64}'))
reads = json.loads(base64.b64decode('${readsB64}'))
originals = [p for p in [unit.get('reads_1'), unit.get('reads_2')] + unit.get('singletons', []) if p]
if len(reads) != len(originals):
    raise ValueError('Staged FASTQ count disagrees with clean manifest')
unit['staged_input_paths'] = dict(zip(originals, reads))
with open('unit.json', 'w') as handle:
    json.dump(unit, handle, sort_keys=True)
PY
    python3 ${code}/reads_assembly.py assemble --unit unit.json --outdir '${id}' --threads ${task.cpus} --memory-gib ${params.assembly_memory_gib}
    """
}
process READS_ANNOTATE {
    tag "$id"
    input:
    tuple val(id), val(material), path(unit_dir)
    path db_manifest
    path code
    val code_digest
    output:
    tuple val(id), val(material), path("${id}_annotated"), emit: units
    script:
    def db = params.bakta_db ? "--db '${params.bakta_db}'" : ''
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    # code_sha256=${code_digest}
    python3 ${code}/reads_assembly.py annotate --unit '${unit_dir}/assembly_unit.json' --outdir '${id}_annotated' --threads ${task.cpus} ${db}
    """
}
process ASSEMBLY_COLLECT {
    publishDir params.outdir, mode: 'copy'
    input:
    path unitdirs
    path code
    val code_digest
    output:
    path 'assemblies', emit: assemblies
    script:
    def dirs = (unitdirs instanceof List ? unitdirs : [unitdirs]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    # code_sha256=${code_digest}
    python3 ${code}/reads_publish.py --assembly-dirs ${dirs} --outdir assemblies --publish-root '${file(params.outdir).toAbsolutePath().normalize()}/assemblies'
    """
}
