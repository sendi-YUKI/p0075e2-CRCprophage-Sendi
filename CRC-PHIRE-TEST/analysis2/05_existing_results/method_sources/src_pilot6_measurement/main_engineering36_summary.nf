nextflow.enable.dsl=2
include { CORE_RESEARCH_SUMMARY } from './modules/local/core_research'
workflow {
    def f=new groovy.json.JsonSlurper().parseText(file(params.engineering_fixture,checkIfExists:true).text)
    if(f.purpose!='engineering_smoke' || f.campaign!='basic_build_20260930' || f.input_kind!='synthetic')error 'Registered artificial campaign only'
    f.files.each { r ->
        def p=file(r.path,checkIfExists:true)
        if(!p.toAbsolutePath().normalize().toString().startsWith('/srv/CRC-PHIRE/projects/crc-pipeline/results/methods_prepare_20261001/'))error 'Fixture outside scope'
        if(java.security.MessageDigest.getInstance('SHA-256').digest(p.bytes).encodeHex().toString()!=r.sha256)error 'Fixture changed'
    }
    def a=file(f.sequence_fixture)
    def helper={n -> file("${projectDir}/bin/${n}",checkIfExists:true)}
    CORE_RESEARCH_SUMMARY(a.resolve('catalog'),a.resolve('locus'),a.resolve('qc'),[a.resolve('genes')],[a.resolve('functions')],a.resolve('function_summary'),a.resolve('search.json'),
        file(f.whole_evidence),true,helper('core_research_summary.py'),helper('catalog_policy.py'),helper('mainline_common.py'),file("${projectDir}/assets/aimbeta/v1/contract.schema.json"),helper('genome_callability.py'))
}
workflow.onComplete {
    def d=new File("${params.outdir}/pipeline_info");d.mkdirs()
    new File(d,'engineering_status.json').text=groovy.json.JsonOutput.toJson([success:workflow.success,exit_status:workflow.exitStatus,session_id:workflow.sessionId.toString(),scope:'Affected summary input only; reused hash-bound synthetic search evidence; no analysis resume'])
}
