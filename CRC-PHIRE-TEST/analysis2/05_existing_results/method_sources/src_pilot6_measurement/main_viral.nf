nextflow.enable.dsl=2
include { VIRAL_DISCOVERY } from './workflows/viral_discovery'

workflow {
    def evidenceTools = (params.viral_evidence_tools ?: '').tokenize(',')
    if (evidenceTools.any { !(it in ['bacphlip','vibrant','vcontact3']) } || evidenceTools.unique(false).size() != evidenceTools.size()) error 'Invalid/duplicate viral_evidence_tools'
    if (evidenceTools && System.getenv('CRC_PHAGE_RUNTIME') != 'spark-native') error 'Viral evidence currently requires prepared Spark native environments'

    if (!params.assembly_manifest) error 'Required: --assembly_manifest assembly_manifest.json'
    if (!(params.catalog_scope in ['host_resolved_prophage','virome_discovery','combined_exploratory'])) error 'Explicit --catalog_scope required'
    if (!params.viral_gene_container || !(params.viral_gene_container ==~ /(?:.+@)?sha256:[a-f0-9]{64}/)) error 'Immutable SHA256 --viral_gene_container is required; dependencies not ready'
    def manifest = file(params.assembly_manifest, checkIfExists:true)
    def doc = new groovy.json.JsonSlurper().parseText(manifest.text)
    if (!(doc.schema_version in [1,'viral_catalog_v1'])) error 'Unsupported or missing assembly manifest schema_version'
    if (!(doc.units instanceof List)) error 'Assembly manifest requires units list (empty list explicitly supported)'
    def ids = []
    doc.units.each { unit ->
        unit.assembly_id = unit.assembly_id ?: unit.unit_id
        unit.unit_id = unit.assembly_id
        if (!(unit.assembly_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]{0,127}/)) error 'Invalid assembly_id/unit_id'
        if (!(unit.sha256 ==~ /[a-f0-9]{64}/)) error "Missing SHA256 for ${unit.assembly_id}"
        if (ids.contains(unit.assembly_id)) error 'Duplicate assembly_id'
        ids << unit.assembly_id
        def fasta = java.nio.file.Paths.get(unit.assembly_fasta)
        unit.assembly_fasta = (fasta.isAbsolute() ? fasta : manifest.parent.resolve(fasta)).normalize().toString()
        if (unit.contig_map) {
            def mp = java.nio.file.Paths.get(unit.contig_map)
            unit.contig_map = (mp.isAbsolute() ? mp : manifest.parent.resolve(mp)).normalize().toString()
        }
    }
    // Validated primary outputs are imported explicitly; callers are not repeated.
    def orderedUnits = doc.units.sort { a,b -> a.assembly_id <=> b.assembly_id }
    def importedFinals = orderedUnits.collect { u -> file(u.primary_final_dir, checkIfExists:true) }
    def importedPrepared = orderedUnits.collect { u -> file(u.primary_prepared_dir, checkIfExists:true) }
    def encoded = groovy.json.JsonOutput.toJson(doc).bytes.encodeBase64().toString()
    def emptyManifest = file("${projectDir}/assets/viral_empty_manifest.json")
    def legacyFastas = []
    if (params.legacy_source_manifest) {
        def lmf = file(params.legacy_source_manifest, checkIfExists:true)
        def ld = new groovy.json.JsonSlurper().parseText(lmf.text)
        legacyFastas = ld.units.sort { a,b -> (a.assembly_id ?: a.unit_id) <=> (b.assembly_id ?: b.unit_id) }.collect { u ->
            def fp = java.nio.file.Paths.get(u.assembly_fasta)
            file(fp.isAbsolute() ? fp : lmf.parent.resolve(fp), checkIfExists:true)
        }
    }
    def researchFastas = []
    if (params.research_catalog) {
        if (!params.research_source_manifest) error 'Research bridge requires research_source_manifest'
        def rm = file(params.research_source_manifest,checkIfExists:true)
        def rd = new groovy.json.JsonSlurper().parseText(rm.text)
        researchFastas = rd.units.sort { a,b -> (a.assembly_id ?: a.unit_id) <=> (b.assembly_id ?: b.unit_id) }.collect { u ->
            def fp = java.nio.file.Paths.get(u.assembly_fasta)
            file(fp.isAbsolute() ? fp : rm.parent.resolve(fp),checkIfExists:true)
        }
    }
    VIRAL_DISCOVERY(Channel.value(importedFinals), Channel.value(importedPrepared), encoded,
        file(params.genomad_db, checkIfExists:true), file(params.checkv_db, checkIfExists:true),
        file(params.functions_db, checkIfExists:true),
        params.legacy_catalog ? file(params.legacy_catalog, checkIfExists:true) : emptyManifest, params.legacy_catalog as boolean,
        params.legacy_source_manifest ? file(params.legacy_source_manifest, checkIfExists:true) : emptyManifest, params.legacy_source_manifest as boolean, legacyFastas,
        params.canonical_gene_tables ? params.canonical_gene_tables.tokenize(',').collect { file(it, checkIfExists:true) } : [],
        params.research_catalog ? file(params.research_catalog,checkIfExists:true) : emptyManifest,
        params.research_source_manifest ? file(params.research_source_manifest,checkIfExists:true) : emptyManifest,
        researchFastas)
}
