nextflow.enable.dsl=2
include { READS } from './workflows/reads'
workflow {
    if (!params.raw_reads_manifest || !params.outdir) error 'Require --raw_reads_manifest and --outdir (bin/run-reads prepares R0)'
    def manifest=file(params.raw_reads_manifest,checkIfExists:true)
    def data=new groovy.json.JsonSlurper().parse(manifest.toFile())
    if(data.schema_version!=1 || !data.units) error 'Invalid or empty raw reads manifest'
    def units=data.units.collect { u -> tuple(u.unit_id,file("${manifest.parent}/units/${u.unit_id}.json",checkIfExists:true),([u.reads_1,u.reads_2]+(u.source_singletons ?: [])).findAll{it}.collect { file(it,checkIfExists:true) }) }
    def codeNames = ['reads_io.py','reads_process.py','reads_assembly.py','reads_publish.py']
    def codeContent = codeNames.collect { n -> n + ':' + java.security.MessageDigest.getInstance('SHA-256').digest(file("${projectDir}/bin/${n}",checkIfExists:true).bytes).encodeHex().toString() }.join('\n')
    def codeDigest = java.security.MessageDigest.getInstance('SHA-256').digest(codeContent.bytes).encodeHex().toString()
    READS(Channel.fromList(units),file(params.preprocessing_profile,checkIfExists:true),manifest,file(params.database_manifest,checkIfExists:true),file(params.software_manifest ?: "${projectDir}/assets/reads_tools.json",checkIfExists:true),file("${projectDir}/bin"),codeDigest)
}
