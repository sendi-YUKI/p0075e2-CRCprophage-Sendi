include { READS_CLEAN; READS_COLLECT } from '../modules/local/reads'
workflow READS {
    take: raw_units; profile; manifest; databases; software; code; code_digest
    main:
    READS_CLEAN(raw_units,profile,databases,software,code,code_digest)
    READS_COLLECT(READS_CLEAN.out.units.collect(),profile,manifest,code,code_digest)
    emit: clean = READS_COLLECT.out.clean
}
