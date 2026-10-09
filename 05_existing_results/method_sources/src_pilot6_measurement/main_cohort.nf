nextflow.enable.dsl=2
include { CRC_PHAGE_COHORT } from './workflows/cohort'
workflow {
    if (!params.cohort_prepared) error 'Run cohort_cli.py prepare first; required --cohort_prepared'
    if (!params.cohort_mapping_container) error 'Fixed B1 mapping container missing; use -c conf/cohort.config and installed assets/reads_tools.json'
    if (!(params.cohort_mapping_container.startsWith('sha256:') || params.cohort_mapping_container.contains('@sha256:')))
        error 'B1 mapping container must be pinned by digest'
    if (params.cohort_enable_metaphlan) {
        if (params.cohort_task_memory_budget_gb < params.cohort_metaphlan_required_memory_gb)
            error "MetaPhlAn requires ${params.cohort_metaphlan_required_memory_gb} GB; declared task budget ${params.cohort_task_memory_budget_gb} GB. Missing ${params.cohort_metaphlan_required_memory_gb - params.cohort_task_memory_budget_gb} GB. No smaller database is substituted."
        if (!params.cohort_metaphlan_container || !params.cohort_metaphlan_db || !params.cohort_metaphlan_index || !params.cohort_metaphlan_database_manifest)
            error 'Enabled MetaPhlAn requires fixed container, full DB/index and database manifest'
    }
    prepared = file(params.cohort_prepared, checkIfExists:true)
    manifest = new groovy.json.JsonSlurper().parse(prepared.resolve('clean_reads_manifest.json').toFile())
    if (manifest.schema_version != 1 || !manifest.containsKey('units')) error 'Invalid clean reads manifest'
    seen = [] as Set
    rows = manifest.units.collect { row ->
        if (!(row.unit_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) || seen.contains(row.unit_id)) error 'Invalid/duplicate unit ID'
        seen.add(row.unit_id)
        unitJson = prepared.resolve("units/${row.unit_id}.json")
        if (new groovy.json.JsonSlurper().parse(unitJson.toFile()) != row) error 'Prepared unit JSON differs from clean reads manifest'
        if (!(row.layout in ['SE','PE']) || row.platform != 'ILLUMINA' || row.molecule != 'DNA') error 'Unsupported readset'
        r2 = row.layout == 'PE' ? [file(row.reads_2, checkIfExists:true)] : []
        tuple([id:row.unit_id, material_type:row.material_type, layout:row.layout],
              prepared.resolve("units/${row.unit_id}.json"),
              file(row.reads_1, checkIfExists:true), r2)
    }
    helpers = Channel.value(file("${projectDir}/bin/cohort_*.py").sort { it.name })
    CRC_PHAGE_COHORT(Channel.fromList(rows), prepared, helpers)
}
workflow.onComplete {
    def directory = new File("${params.outdir}/pipeline_info")
    directory.mkdirs()
    new File(directory, 'workflow_status.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        success:workflow.success, exit_status:workflow.exitStatus, session_id:workflow.sessionId.toString(),
        current_target:'READS_B1_B2_BUILD', workflow:'B1', real_data_validated:false, scientifically_calibrated:false
    ])) + '\n'
}
