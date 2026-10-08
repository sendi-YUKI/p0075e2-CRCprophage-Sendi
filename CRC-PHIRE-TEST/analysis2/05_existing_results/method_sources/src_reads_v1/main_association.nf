nextflow.enable.dsl=2
include { ASSOCIATION } from './workflows/association'

workflow {
    if (!params.association_long || !params.association_metadata || !params.association_spec)
        error 'Required --association_long --association_metadata --association_spec'
    if (!params.association_container)
        error 'Association container unavailable: run bootstrap/association_container.sh'
    ASSOCIATION(
        file(params.association_long, checkIfExists:true),
        file(params.association_metadata, checkIfExists:true),
        file(params.association_spec, checkIfExists:true))
}

workflow.onComplete {
    def dir = new File("${params.outdir}/pipeline_info")
    dir.mkdirs()
    new File(dir,'workflow_status.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        success:workflow.success, exit_status:workflow.exitStatus, session_id:workflow.sessionId.toString(),
        run_name:workflow.runName, current_target:'READS_B1_B2_BUILD', workflow_entry:'association',
        completed:workflow.complete.toString(), real_data_validated:false, scientifically_calibrated:false
    ])) + '\n'
}
