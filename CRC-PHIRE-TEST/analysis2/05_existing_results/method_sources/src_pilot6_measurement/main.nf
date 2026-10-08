#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { CRC_PHAGE } from './workflows/crc_phage'

workflow {
    // Native execution must come through the existing allocation/lock/prefix-checked
    // spark-run chain. Plain host PATH remains unsupported; task-shell keeps the guard.
    def nativeRuntime = System.getenv('CRC_PHAGE_RUNTIME') == 'spark-native' &&
        workflow.profile.tokenize(',').contains('spark_native') &&
        System.getenv('SLURM_JOB_ID') && System.getenv('CRC_SPARK_HEAVY_LOCK_OWNER')
    if (workflow.containerEngine != 'docker' && !nativeRuntime)
        error 'A0+A1 requires Docker or the allocation-checked spark-run native runtime; plain host-PATH execution is disabled.'
    if (!params.input) error 'Required: --input /path/to/genomes.tsv'
    if (!params.genomad_db || !params.checkv_db) error 'Required: --genomad_db and --checkv_db'
    def sheet = file(params.input, checkIfExists: true)
    if (!sheet.isFile()) error "Samplesheet is not a regular file: ${sheet}"
    def table = sheet.splitCsv(header:false, sep:'\t', strip:false)
    if (!table || table.size()<2) error 'Input samplesheet has no genomes'
    def headers = table[0].collect { it ?: '' }
    headers[0] = headers[0].replaceFirst('^\uFEFF','')
    if (headers.any{ !it } || headers.toSet().size()!=headers.size()) error 'Empty or duplicate samplesheet header'
    if (!headers.containsAll(['genome_id','fasta'])) error 'Input columns must include genome_id and fasta'
    def rows = table.drop(1).collect { values ->
        if (values.size()!=headers.size() || values.any { it==null }) error 'Samplesheet row has the wrong number of columns'
        headers.withIndex().collectEntries { name, index -> [(name):values[index]] }
    }
    def ids = [] as Set
    def source_paths = [] as Set
    def samples = rows.collect { row ->
        if (!row.containsKey('genome_id') || !row.containsKey('fasta')) error 'Input columns must include genome_id and fasta'
        def id = row.genome_id
        if (!id || !(id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]{0,127}/)) error "Invalid genome_id: ${id}"
        if (!ids.add(id)) error "Duplicate genome_id: ${id}"
        if (!row.fasta || row.fasta!=row.fasta.trim()) error "Missing or whitespace-padded fasta for ${id}"
        def input_name = row.fasta.startsWith('~/') ? System.getenv('HOME') + row.fasta.substring(1) : row.fasta
        def source = java.nio.file.Paths.get(input_name)
        def resolved = source.isAbsolute() ? source : sheet.parent.resolve(source)
        def fasta = file(resolved.toString(), checkIfExists:true)
        if (!fasta.isFile() || fasta.size()==0) error "FASTA is not a nonempty regular file for ${id}: ${fasta}"
        if (!source_paths.add(fasta.toRealPath().toString())) error "Same FASTA assigned to multiple genome IDs: ${fasta}"
        tuple([id:id], fasta)
    }
    def genomad_db = file(params.genomad_db, checkIfExists:true)
    def checkv_db = file(params.checkv_db, checkIfExists:true)
    if (!genomad_db.isDirectory() || !checkv_db.isDirectory()) error 'Database inputs must be directories.'
    CRC_PHAGE(Channel.fromList(samples), genomad_db, checkv_db)
}

workflow.onComplete {
    new File("${params.outdir}/pipeline_info").mkdirs()
    new File("${params.outdir}/pipeline_info/workflow_status.json").text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        success:workflow.success, exit_status:workflow.exitStatus, error_message:workflow.errorMessage,
        run_name:workflow.runName, session_id:workflow.sessionId.toString(), nextflow_version:nextflow.version.toString(),
        command_line:workflow.commandLine, start:workflow.start.toString(), complete:workflow.complete.toString(),
        work_dir:workflow.workDir.toString(), parameters:params
    ]))
}
