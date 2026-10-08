nextflow.enable.dsl=2
include { RESEARCH_V1_BRIDGE } from './modules/local/research_bridge'
include { CORE_WHOLE_GENOME_EVIDENCE } from './modules/local/core_research'
include { HOST_SOURCE_ASSOCIATION } from './modules/local/host_source'
include { HOST_RECOMBINATION } from './modules/local/host_recombination'
include { VIRAL_KOFAM } from './modules/local/viral_kofam'
include { CRC_PHAGE_COHORT } from './workflows/cohort'
workflow {
    if (!params.engineering_fixture) error 'Explicit campaign fixture registry required'
    def mf=file(params.engineering_fixture,checkIfExists:true)
    def f=new groovy.json.JsonSlurper().parseText(mf.text)
    if(f.purpose!='engineering_smoke' || f.campaign!='basic_build_20260930' || f.input_kind!='synthetic') error 'Only registered synthetic fixtures supported'
    f.files.each { r ->
        def p=file(r.path,checkIfExists:true)
        def actual=java.security.MessageDigest.getInstance('SHA-256').digest(p.bytes).encodeHex().toString()
        if(actual!=r.sha256)error "Fixture changed: ${r.path}"
    }
    if (!params.recombination_only) {
    def seq=file(f.sequence_fixture,checkIfExists:true)
    def bh=['research_bridge.py','mainline_common.py','viral_catalog.py','core_votu.py','phageflow.py'].collect { file("${projectDir}/bin/${it}") }
    RESEARCH_V1_BRIDGE(seq.resolve('research'),seq.resolve('independent'),seq.resolve('sources.json'),[seq.resolve('prepared_G1/genome.fna')],bh)
    def wh=['whole_genome_evidence.py','core_locus.py','catalog_policy.py','core_votu.py','mainline_common.py','phageflow.py'].collect { file("${projectDir}/bin/${it}") }
    CORE_WHOLE_GENOME_EVIDENCE([seq.resolve('prepared_G1')],seq.resolve('research'),seq.resolve('search.json'),wh)
    def hs=file(f.host_fixture,checkIfExists:true)
    HOST_SOURCE_ASSOCIATION(hs.resolve('fixture/research'),hs.resolve('fixture/mask'),hs.resolve('fixture/tree'),
       hs.resolve('fixture/config.json'),hs.resolve('fixture/bundle'),file("${projectDir}/bin/host_source.py"),file("${projectDir}/bin/host_source_models.R"),[],
       file("${projectDir}/bin/host_sensitivity.R"))
    }
    def rc=file(f.recombination_fixture,checkIfExists:true)
    HOST_RECOMBINATION(rc.resolve('config.json'),rc.resolve('assemblies.json'),(0..3).collect { rc.resolve("G${it}.fna") },file("${projectDir}/bin/host_recombination.py"))
    if (!params.recombination_only) {
    def seq=file(f.sequence_fixture,checkIfExists:true)
    VIRAL_KOFAM(file(f.empty_genes),seq.resolve('bridge'),file(f.kofam_database),file(f.kofam_receipt),
      ['core_kofam.py','viral_kofam.py','mainline_common.py','kofam_hmmsearch_serial'].collect { file("${projectDir}/bin/${it}") })
    def prepared=seq.resolve('cohort_prepared')
    def unit=new groovy.json.JsonSlurper().parseText(prepared.resolve('units/U1.json').text)
    CRC_PHAGE_COHORT(Channel.of(tuple([id:'U1',material_type:unit.material_type,layout:unit.layout],prepared.resolve('units/U1.json'),file(unit.reads_1),[])),
       prepared,Channel.value(file("${projectDir}/bin/cohort_*.py").sort { it.name }))
    }
}
workflow.onComplete {
    def folder=new File("${params.outdir}/pipeline_info");folder.mkdirs()
    new File(folder,'engineering_status.json').text=groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
      purpose:'engineering_smoke',success:workflow.success,exit_status:workflow.exitStatus,session_id:workflow.sessionId.toString(),
      biological_input:false,scientific_calibration:false,skipped:'geNomad/CheckV/Panaroo/annotation biological searches; KOfam empty-ORF interface only'
    ]))
}
