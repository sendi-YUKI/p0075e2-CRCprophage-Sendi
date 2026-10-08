include { ASSOCIATION_MODELS } from '../modules/local/association_models'

workflow ASSOCIATION {
    take:
    long_table
    metadata
    spec
    main:
    ASSOCIATION_MODELS(long_table, metadata, spec,
        file("${projectDir}/assets/association_spec.schema.json"),
        file("${projectDir}/bin/association.py"),
        file("${projectDir}/bin/association_models.R"))
    emit:
    result = ASSOCIATION_MODELS.out.result
}
