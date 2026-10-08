include { HOST_RECOMBINATION } from '../modules/local/host_recombination'
include { CORE_GENES; CORE_CHECKM2; CORE_FASTANI; CORE_TAXONOMY_NORMALIZE } from '../modules/local/core_genomes'
include { CORE_QC_COLLECT; CORE_BLAST; CORE_VOTU } from '../modules/local/core_votu'
include { CORE_FUNCTIONS; CORE_FUNCTION_SUMMARY; CORE_ELEMENT_COLLECT } from '../modules/local/core_functions'
include { CORE_PHISPY; CORE_PHAGEBOOST; CORE_PHISPY_UNASSESSED; CORE_COMPARE_CALLERS } from '../modules/local/core_callers'
include { CORE_REPORT } from '../modules/local/core_report'
include { CORE_KOFAM_DB_VALIDATE; CORE_KOFAM; CORE_KOFAM_SUMMARY } from '../modules/local/core_kofam'
include { HOST_CORE_GFF; HOST_CORE_PANAROO; HOST_CORE_ALIGN; HOST_TREE_PRELIMINARY; HOST_MOBILE_MASK; HOST_TREE_FINAL } from '../modules/local/host_phylogeny'

include { CORE_WHOLE_GENOME_EVIDENCE; CORE_RESEARCH_CATALOG; CORE_LOCUS_EVIDENCE; CORE_RESEARCH_SUMMARY } from '../modules/local/core_research'

include { HOST_SOURCE_ASSOCIATION } from '../modules/local/host_source'

workflow CRC_PHAGE_CORE {
    take:
    prepared
    a1
    genbanks
    checkm2_db
    functions_db
    reference_dir
    metadata

    main:
    genomeHelper = file("${projectDir}/bin/core_genomes.py",checkIfExists:true)
    phageHelper = file("${projectDir}/bin/phageflow.py",checkIfExists:true)
    votuHelper = file("${projectDir}/bin/core_votu.py",checkIfExists:true)
    functionsHelper = file("${projectDir}/bin/core_functions.py",checkIfExists:true)
    callerHelper = file("${projectDir}/bin/core_callers.py",checkIfExists:true)
    reportHelper = file("${projectDir}/bin/core_report.py",checkIfExists:true)
    master = a1.resolve('prophage_master.tsv')
    fasta = a1.resolve('candidates.fna')
    qcPolicy = [completeness_min:params.completeness_min, contamination_max:params.contamination_max,
        high_completeness_min:params.high_completeness_min, high_contamination_max:params.high_contamination_max,
        fragmentation_n50:params.fragmentation_n50]
    taxPolicy = [ani_min:params.taxonomy_ani_min, fragment_fraction_min:params.taxonomy_fragment_fraction_min]
    votuPolicy = [ani_percent:params.votu_ani_percent, af_shorter_percent:params.votu_af_shorter_percent]
    CORE_GENES(prepared, genomeHelper, phageHelper)
    CORE_CHECKM2(prepared, checkm2_db, metadata, qcPolicy, genomeHelper, phageHelper)
    CORE_FASTANI(prepared, reference_dir)
    CORE_TAXONOMY_NORMALIZE(CORE_FASTANI.out.raw, reference_dir, metadata, taxPolicy, genomeHelper, phageHelper)
    qcDirs = CORE_CHECKM2.out.qc.collect(flat:false).map { items -> items.sort { a,b -> a[0].id <=> b[0].id }.collect { it[1] } }
    CORE_QC_COLLECT(qcDirs, genomeHelper, phageHelper)
    CORE_BLAST(fasta)
    mainlineCommon = file("${projectDir}/bin/mainline_common.py", checkIfExists:true)
    CORE_FUNCTIONS(CORE_GENES.out.genes, master, functions_db, functionsHelper)
    functionDirs = CORE_FUNCTIONS.out.functions.collect(flat:false).map { items -> items.sort { a,b -> a[0].id <=> b[0].id }.collect { it[1] } }
    CORE_ELEMENT_COLLECT(master, functionDirs, functionsHelper)
    elementClass = CORE_ELEMENT_COLLECT.out.elements.map { it.resolve('element_class.tsv') }
    coreCatalog = Channel.empty()
    locusOutputs = Channel.value([])
    researchOutputs = Channel.value([])
    if (params.catalog_policy) {
        catalogPolicy = file(params.catalog_policy, checkIfExists:true)
        catalogEvidence = file(params.catalog_evidence, checkIfExists:true)
        catalogHelper = file("${projectDir}/bin/catalog_policy.py",checkIfExists:true)
        CORE_RESEARCH_CATALOG(master, fasta, CORE_BLAST.out.blast, CORE_QC_COLLECT.out.qc, elementClass, votuHelper, catalogHelper, mainlineCommon, catalogPolicy, catalogEvidence)
        coreCatalog = CORE_RESEARCH_CATALOG.out.catalog
        preparedDirs = prepared.collect(flat:false).map { items -> items.sort { a,b -> a[0].id <=> b[0].id }.collect { it[1] } }
        CORE_LOCUS_EVIDENCE(preparedDirs, coreCatalog, catalogPolicy, catalogEvidence, file("${projectDir}/bin/core_locus.py",checkIfExists:true), catalogHelper, votuHelper, mainlineCommon)
        locusOutputs = CORE_LOCUS_EVIDENCE.out.locus
        if (params.whole_genome_search_config) {
            CORE_WHOLE_GENOME_EVIDENCE(preparedDirs,coreCatalog,file(params.whole_genome_search_config,checkIfExists:true),
                ['whole_genome_evidence.py','genome_callability.py','core_locus.py','catalog_policy.py','core_votu.py','mainline_common.py','phageflow.py'].collect { file("${projectDir}/bin/${it}",checkIfExists:true) })
        }

    } else {
        CORE_VOTU(master, fasta, CORE_BLAST.out.blast, CORE_QC_COLLECT.out.qc, elementClass, votuHelper, votuPolicy)
        coreCatalog = CORE_VOTU.out.catalog
    }

    CORE_FUNCTION_SUMMARY(coreCatalog.map { it.resolve('vOTU_members.tsv') }, functionDirs, functionsHelper)
    phiInputs = prepared.branch { meta, dir ->
        assessed: genbanks.containsKey(meta.id)
        unassessed: true
    }
    CORE_PHISPY(phiInputs.assessed.map { meta,dir -> tuple(meta,dir,genbanks[meta.id]) }, callerHelper)
    CORE_PHISPY_UNASSESSED(phiInputs.unassessed, callerHelper)
    CORE_PHAGEBOOST(prepared, callerHelper)
    callerDirs = CORE_PHISPY.out.result.mix(CORE_PHAGEBOOST.out.result).mix(CORE_PHISPY_UNASSESSED.out.result)
        .collect(flat:false).map { items -> items.sort { a,b -> "${a[0].id}:${a[1].name}" <=> "${b[0].id}:${b[1].name}" }.collect { it[1] } }
    CORE_COMPARE_CALLERS(master, callerDirs, callerHelper)
    geneDirs = CORE_GENES.out.genes.collect(flat:false).map { items -> items.sort { a,b -> a[0].id <=> b[0].id }.collect { it[1] } }
    taxonomyDirs = CORE_TAXONOMY_NORMALIZE.out.taxonomy.collect(flat:false).map { items -> items.sort { a,b -> a[0].id <=> b[0].id }.collect { it[1] } }
    kofamSummary = Channel.value([])
    hostOutputs = Channel.value([])
    hostSourceOutputs = Channel.value([])
    recombinationOutputs = Channel.value([])
    if (params.enable_kofam) {
        kofamHelper = file("${projectDir}/bin/core_kofam.py", checkIfExists:true)
        kofamShim = file("${projectDir}/bin/kofam_hmmsearch_serial", checkIfExists:true)
        kofamDb = file(params.kofam_db, checkIfExists:true)
        kofamIdentity = new groovy.json.JsonSlurper().parseText(file("${projectDir}/assets/spark_databases.json").text).databases.kofam.manifest_sha256
        CORE_KOFAM_DB_VALIDATE(kofamDb, kofamIdentity, kofamHelper, mainlineCommon)
        CORE_KOFAM(CORE_GENES.out.genes, kofamDb, CORE_KOFAM_DB_VALIDATE.out.receipt, kofamHelper, mainlineCommon, kofamShim)
        kofamDirs = CORE_KOFAM.out.result.collect(flat:false).map { items -> items.sort { a,b -> a[0].id <=> b[0].id }.collect { it[1] } }
        CORE_KOFAM_SUMMARY(coreCatalog.map { it.resolve('vOTU_members.tsv') }, kofamDirs, functionDirs, kofamHelper, mainlineCommon)
        kofamSummary = CORE_KOFAM_SUMMARY.out.summary
    }
    if (params.host_config) {
        hostConfig = file(params.host_config, checkIfExists:true)
        hostPlan = new groovy.json.JsonSlurper().parseText(hostConfig.text)
        hostHelper = file("${projectDir}/bin/host_phylogeny.py", checkIfExists:true)
        maskBundle = file(params.host_mask_bundle, checkIfExists:true)
        hostInput = CORE_GENES.out.genes.join(prepared).filter { meta,genes,genome -> hostPlan.genome_ids.contains(meta.id) }
        HOST_CORE_GFF(hostInput, hostConfig, hostHelper, mainlineCommon)
        hostDirs = HOST_CORE_GFF.out.adapted.collect(flat:false).map { items -> items.sort { a,b -> a[0].id <=> b[0].id }.collect { it[1] } }
        HOST_CORE_PANAROO(hostDirs, CORE_QC_COLLECT.out.qc.map { it.resolve('genome_qc.tsv') }, hostConfig, hostHelper, mainlineCommon)
        HOST_CORE_ALIGN(HOST_CORE_PANAROO.out.clusters, hostConfig, hostHelper, mainlineCommon)
        HOST_TREE_PRELIMINARY(HOST_CORE_ALIGN.out.alignment, hostConfig, hostHelper, mainlineCommon)
        HOST_MOBILE_MASK(HOST_CORE_ALIGN.out.alignment, HOST_TREE_PRELIMINARY.out.tree, master, maskBundle, hostConfig, hostHelper, mainlineCommon)
        HOST_TREE_FINAL(HOST_MOBILE_MASK.out.masked, hostConfig, hostHelper, mainlineCommon)
        hostOutputs = HOST_CORE_ALIGN.out.alignment.combine(HOST_TREE_PRELIMINARY.out.tree).combine(HOST_MOBILE_MASK.out.masked).combine(HOST_TREE_FINAL.out.tree)
    }
    if (params.catalog_policy) {
        CORE_RESEARCH_SUMMARY(coreCatalog, locusOutputs, CORE_QC_COLLECT.out.qc, geneDirs, functionDirs, CORE_FUNCTION_SUMMARY.out.summary, kofamSummary,
            params.whole_genome_search_config ? CORE_WHOLE_GENOME_EVIDENCE.out.evidence : file("${projectDir}/assets/whole_genome_search.proposed.json"),
            params.whole_genome_search_config as boolean,
            file("${projectDir}/bin/core_research_summary.py",checkIfExists:true), catalogHelper, mainlineCommon, file("${projectDir}/assets/aimbeta/v1/contract.schema.json",checkIfExists:true), file("${projectDir}/bin/genome_callability.py",checkIfExists:true))
        researchOutputs = CORE_RESEARCH_SUMMARY.out.summary
    }
    if (params.host_recombination_config) {
        rcfg=file(params.host_recombination_config,checkIfExists:true)
        rplan=new groovy.json.JsonSlurper().parseText(rcfg.text)
        rmanifest=file(rplan.assembly_manifest,checkIfExists:true)
        rsource=new groovy.json.JsonSlurper().parseText(rmanifest.text)
        rfastas=rsource.units.sort { a,b -> (a.assembly_id ?: a.unit_id) <=> (b.assembly_id ?: b.unit_id) }.collect { u ->
            def pp=java.nio.file.Paths.get(u.assembly_fasta)
            file(pp.isAbsolute() ? pp : rmanifest.parent.resolve(pp),checkIfExists:true)
        }
        HOST_RECOMBINATION(rcfg,rmanifest,rfastas,file("${projectDir}/bin/host_recombination.py"))
        recombinationOutputs=HOST_RECOMBINATION.out.result
    }
    if (params.host_association_config) {
        if (!params.host_config || !params.catalog_policy) error 'Host-source association requires host and research catalog branches'
        HOST_SOURCE_ASSOCIATION(researchOutputs, params.host_recombination_config ? HOST_RECOMBINATION.out.masked : HOST_MOBILE_MASK.out.masked, params.host_recombination_config ? HOST_RECOMBINATION.out.tree : HOST_TREE_FINAL.out.tree,
            file(params.host_association_config,checkIfExists:true), file(params.host_association_bundle,checkIfExists:true),
            file("${projectDir}/bin/host_source.py",checkIfExists:true), file("${projectDir}/bin/host_source_models.R",checkIfExists:true), recombinationOutputs, file("${projectDir}/bin/host_sensitivity.R",checkIfExists:true))
        hostSourceOutputs = HOST_SOURCE_ASSOCIATION.out.result
    }
    CORE_REPORT(a1, coreCatalog, CORE_QC_COLLECT.out.qc, geneDirs, taxonomyDirs,
        functionDirs, CORE_FUNCTION_SUMMARY.out.summary, CORE_COMPARE_CALLERS.out.result, kofamSummary, hostOutputs, locusOutputs, researchOutputs, hostSourceOutputs, reportHelper)

    emit:
    report = CORE_REPORT.out.report
}
