include { COHORT_INDEX; COHORT_MAPPABILITY; COHORT_MEASURE; COHORT_METAPHLAN; COHORT_REPORT } from '../modules/local/cohort'
workflow CRC_PHAGE_COHORT {
    take:
    units
    prepared
    helpers
    main:
    profile = prepared.resolve('measurement_profile.json')
    detection = prepared.resolve('detection_profiles.json')
    COHORT_INDEX(prepared.resolve('reference'), helpers)
    COHORT_MAPPABILITY(COHORT_INDEX.out.index, profile, helpers)
    COHORT_MEASURE(units, COHORT_INDEX.out.index, profile, detection, COHORT_MAPPABILITY.out.summary, helpers)
    hostResults = Channel.empty()
    if (params.cohort_enable_metaphlan) {
        COHORT_METAPHLAN(units.filter { meta, unit, read1, read2 -> meta.material_type == 'bulk_metagenome' },
                        file(params.cohort_metaphlan_db, checkIfExists:true),
                        file(params.cohort_metaphlan_database_manifest, checkIfExists:true),
                        params.cohort_metaphlan_index, helpers)
        hostResults = COHORT_METAPHLAN.out.results
    }
    COHORT_REPORT(prepared, COHORT_MEASURE.out.results.map { meta, result -> result }.collect().ifEmpty([]),
                  hostResults.collect().ifEmpty([]), helpers, params.cohort_metaphlan_disabled_reason ?: 'not_enabled')
    emit:
    report = COHORT_REPORT.out.report
}
