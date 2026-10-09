include { CORE_KOFAM_DB_VALIDATE as VIRAL_KOFAM_DB_VALIDATE } from '../modules/local/core_kofam'
include { VIRAL_KOFAM } from '../modules/local/viral_kofam'
include { RESEARCH_V1_BRIDGE } from '../modules/local/research_bridge'
include { GENOMAD_ENDTOEND as VIRAL_GENOMAD } from '../modules/nf-core/genomad/endtoend/main'
include { CHECKV_ENDTOEND as VIRAL_CHECKV } from '../modules/nf-core/checkv/endtoend/main'
include { NORMALIZE_CALLS as VIRAL_NORMALIZE; ANNOTATE_QUALITY as VIRAL_QUALITY; EMPTY_QUALITY as VIRAL_EMPTY_QUALITY } from '../modules/local/a0_a1'
include { VIRAL_PREPARE; VIRAL_COLLECT; VIRAL_BLAST; VIRAL_CLUSTER; VIRAL_GENES; VIRAL_GENE_CROSSWALK; VIRAL_FUNCTIONS; VIRAL_FUNCTION_SUMMARY } from '../modules/local/viral'

include { VIRAL_BACPHLIP; VIRAL_VIBRANT; VIRAL_VCONTACT3; VIRAL_EVIDENCE_SUMMARY } from '../modules/local/viral_evidence'

workflow VIRAL_DISCOVERY {
    take:
    finals
    prepared
    manifest_b64
    genomad_db
    checkv_db
    functions_db
    legacy_catalog
    use_legacy
    legacy_sources
    use_legacy_sources
    legacy_fastas
    canonical_genes
    research_catalog
    research_sources
    research_fastas
    main:
    helpers = ['viral_catalog.py','viral_genes.py','viral_functions.py','phageflow.py','core_votu.py','core_functions.py'].collect { file("${projectDir}/bin/${it}", checkIfExists:true) }
    base_helper = file("${projectDir}/bin/phageflow.py", checkIfExists:true)
    VIRAL_COLLECT(finals, prepared, manifest_b64, legacy_catalog, use_legacy, legacy_sources, use_legacy_sources, legacy_fastas, helpers)
    VIRAL_BLAST(VIRAL_COLLECT.out.catalog)
    VIRAL_CLUSTER(VIRAL_COLLECT.out.catalog, VIRAL_BLAST.out.alignments, helpers)
    if (params.research_catalog) {
        bridge_helpers = ['research_bridge.py','mainline_common.py','viral_catalog.py','core_votu.py','phageflow.py'].collect { file("${projectDir}/bin/${it}",checkIfExists:true) }
        RESEARCH_V1_BRIDGE(research_catalog,VIRAL_CLUSTER.out.catalog,research_sources,research_fastas,bridge_helpers)
    }

    VIRAL_GENES(VIRAL_COLLECT.out.catalog, helpers, file("${projectDir}/assets/viral_gene_tools.json", checkIfExists:true))
    if (params.enable_viral_kofam) {
        VIRAL_KOFAM_DB_VALIDATE(file(params.kofam_db,checkIfExists:true),params.kofam_db_manifest_sha256,
            file("${projectDir}/bin/core_kofam.py"),file("${projectDir}/bin/mainline_common.py"))
        VIRAL_KOFAM(VIRAL_GENES.out.genes,VIRAL_CLUSTER.out.catalog,file(params.kofam_db),
            VIRAL_KOFAM_DB_VALIDATE.out.receipt,['core_kofam.py','viral_kofam.py','mainline_common.py','kofam_hmmsearch_serial'].collect { file("${projectDir}/bin/${it}",checkIfExists:true) })
    }
    VIRAL_GENE_CROSSWALK(VIRAL_COLLECT.out.catalog, VIRAL_GENES.out.genes, canonical_genes, helpers)
    VIRAL_FUNCTIONS(VIRAL_COLLECT.out.catalog, VIRAL_GENES.out.genes, functions_db, file("${projectDir}/assets/core_functions_tools.json"), helpers)
    VIRAL_FUNCTION_SUMMARY(VIRAL_CLUSTER.out.catalog, VIRAL_FUNCTIONS.out.functions, helpers)
    evidence_tools = (params.viral_evidence_tools ?: '').tokenize(',')
    evidence_helper = file("${projectDir}/bin/viral_evidence.py", checkIfExists:true)
    evidence_helpers = ['viral_evidence.py','evidence_crosswalk.py','viral_evidence_summary.py'].collect { file("${projectDir}/bin/${it}", checkIfExists:true) }
    evidence_contract = file("${projectDir}/assets/viral_evidence_tools.json", checkIfExists:true)
    evidence_outputs = Channel.empty()
    if ('bacphlip' in evidence_tools) {
        VIRAL_BACPHLIP('bacphlip', VIRAL_CLUSTER.out.catalog, evidence_contract, evidence_helper)
        evidence_outputs = evidence_outputs.mix(VIRAL_BACPHLIP.out.evidence)
    }
    if ('vibrant' in evidence_tools) {
        VIRAL_VIBRANT('vibrant', VIRAL_CLUSTER.out.catalog, evidence_contract, evidence_helpers, VIRAL_GENES.out.genes, canonical_genes)
        evidence_outputs = evidence_outputs.mix(VIRAL_VIBRANT.out.evidence)
    }
    if ('vcontact3' in evidence_tools) {
        VIRAL_VCONTACT3('vcontact3', VIRAL_CLUSTER.out.catalog, evidence_contract, evidence_helper)
        evidence_outputs = evidence_outputs.mix(VIRAL_VCONTACT3.out.evidence)
    }
    VIRAL_EVIDENCE_SUMMARY(VIRAL_CLUSTER.out.catalog, evidence_outputs.collect().ifEmpty([]), evidence_tools.join(','), evidence_helpers)
    emit:
    catalog = VIRAL_CLUSTER.out.catalog
    genes = VIRAL_GENES.out.genes
    gene_crosswalk = VIRAL_GENE_CROSSWALK.out.crosswalk
    functions = VIRAL_FUNCTIONS.out.functions
    function_summary = VIRAL_FUNCTION_SUMMARY.out.summary
    evidence_summary = VIRAL_EVIDENCE_SUMMARY.out.summary
}
