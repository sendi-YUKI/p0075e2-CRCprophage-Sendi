nextflow.enable.dsl=2
include { CORE_WHOLE_GENOME_EVIDENCE; CORE_RESEARCH_SUMMARY } from './modules/local/core_research'
include { HOST_RECOMBINATION } from './modules/local/host_recombination'
include { CRC_PHAGE_COHORT } from './workflows/cohort'
workflow {
    def f=new groovy.json.JsonSlurper().parseText(file(params.engineering_fixture,checkIfExists:true).text)
    if(f.purpose!='engineering_smoke' || f.campaign!='basic_build_20260930' || f.input_kind!='synthetic') error 'Registered artificial campaign only'
    f.files.each { r ->
        def p=file(r.path,checkIfExists:true)
        if(!p.toAbsolutePath().normalize().toString().startsWith('/srv/CRC-PHIRE/projects/crc-pipeline/results/methods_prepare_20261001/'))error 'Fixture outside authorized scope'
        if(java.security.MessageDigest.getInstance('SHA-256').digest(p.bytes).encodeHex().toString()!=r.sha256)error 'Fixture hash changed'
    }
    def a=file(f.sequence_fixture);def rc=file(f.recombination_fixture);def prepared=file(f.prepared)
    params.whole_genome_search_config=a.resolve('search.json').toString()
    params.enable_kofam=false
    def helper={ n -> file("${projectDir}/bin/${n}",checkIfExists:true) }
    CORE_WHOLE_GENOME_EVIDENCE((1..4).collect { a.resolve("prepared_G${it}") },a.resolve('catalog'),a.resolve('search.json'),
        ['whole_genome_evidence.py','genome_callability.py','core_locus.py','catalog_policy.py','core_votu.py','mainline_common.py','phageflow.py'].collect(helper))
    CORE_RESEARCH_SUMMARY(a.resolve('catalog'),a.resolve('locus'),a.resolve('qc'),[a.resolve('genes')],[a.resolve('functions')],a.resolve('function_summary'),a.resolve('search.json'),
        CORE_WHOLE_GENOME_EVIDENCE.out.evidence,true,helper('core_research_summary.py'),helper('catalog_policy.py'),helper('mainline_common.py'),file("${projectDir}/assets/aimbeta/v1/contract.schema.json"),helper('genome_callability.py'))
    HOST_RECOMBINATION(rc.resolve('config.json'),rc.resolve('assemblies.json'),(0..3).collect { rc.resolve("G${it}.fna") },helper('host_recombination.py'))
    def units=['SE','PE'].collect { name ->
        def u=new groovy.json.JsonSlurper().parseText(prepared.resolve("units/${name}.json").text)
        tuple([id:name,material_type:u.material_type,layout:u.layout],prepared.resolve("units/${name}.json"),file(u.reads_1),u.reads_2 ? file(u.reads_2) : [])
    }
    CRC_PHAGE_COHORT(Channel.fromList(units),prepared,Channel.value(file("${projectDir}/bin/cohort_*.py").sort { it.name }))
}
workflow.onComplete {
    def d=new File("${params.outdir}/pipeline_info");d.mkdirs()
    new File(d,'engineering_status.json').text=groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        success:workflow.success,exit_status:workflow.exitStatus,session_id:workflow.sessionId.toString(),purpose:'engineering_smoke',
        biological_inputs:false,scientific_calibration:false,scope:'Changed callability/research summary, recombination and SE/PE reference measurement/report modules only']))
}
