// Repository-root entry preserves projectDir for the actual native profile.
nextflow.enable.dsl=2
include { NATIVE_CACHE_CHECK } from './tests/spark/native_cache/main'
workflow { NATIVE_CACHE_CHECK() }
workflow.onComplete {
    new File(params.cache_probe_status).text = groovy.json.JsonOutput.toJson([
        success:workflow.success,exit_status:workflow.exitStatus,
        uuid:workflow.sessionId.toString(),work_dir:workflow.workDir.toString()
    ]) + '\n'
}
