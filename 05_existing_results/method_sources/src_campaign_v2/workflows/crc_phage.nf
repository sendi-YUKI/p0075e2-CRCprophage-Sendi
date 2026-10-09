include { GENOMAD_ENDTOEND } from '../modules/nf-core/genomad/endtoend/main'
include { CHECKV_ENDTOEND } from '../modules/nf-core/checkv/endtoend/main'
include { PREPARE_GENOME; NORMALIZE_CALLS; ANNOTATE_QUALITY; EMPTY_QUALITY; COLLECT_REPORT } from '../modules/local/a0_a1'

workflow CRC_PHAGE {
    take:
    genomes
    genomad_db
    checkv_db

    main:
    // Explicit script staging makes container access and resume invalidation
    // independent of an external host PATH or an implicit project bind mount.
    helper = file("${projectDir}/bin/phageflow.py", checkIfExists: true)
    PREPARE_GENOME(genomes, helper)
    GENOMAD_ENDTOEND(PREPARE_GENOME.out.prepared.map { meta, dir -> tuple(meta, dir.resolve('genome.fna')) }, tuple([id:'genomad_db_1.9'], genomad_db))
    NORMALIZE_CALLS(PREPARE_GENOME.out.prepared.join(GENOMAD_ENDTOEND.out.genomad_results, failOnDuplicate:true, failOnMismatch:true), helper)
    selected = NORMALIZE_CALLS.out.normalized.branch { meta, dir ->
        positive: dir.resolve('candidate_count.txt').text.trim().toInteger() > 0
        empty: true
    }
    CHECKV_ENDTOEND(selected.positive.map { meta, dir -> tuple(meta,dir.resolve('candidates.fna')) }, checkv_db)
    quality_inputs = selected.positive.join(CHECKV_ENDTOEND.out.quality_summary, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, dir, quality -> tuple(meta,dir,quality.parent) }
    ANNOTATE_QUALITY(quality_inputs, helper)
    EMPTY_QUALITY(selected.empty, helper)
    // Task completion order can differ on resume. Preserve sample metadata
    // until sorting so the collector receives the same ordered path lists.
    final_dirs = ANNOTATE_QUALITY.out.finished.mix(EMPTY_QUALITY.out.finished)
        .collect(flat:false)
        .map { entries -> entries.sort { left, right -> left[0].id <=> right[0].id }.collect { entry -> entry[1] } }
    prepared_dirs = PREPARE_GENOME.out.prepared
        .collect(flat:false)
        .map { entries -> entries.sort { left, right -> left[0].id <=> right[0].id }.collect { entry -> entry[1] } }
    COLLECT_REPORT(final_dirs, prepared_dirs, helper)

    emit:
    report = COLLECT_REPORT.out.report
}
