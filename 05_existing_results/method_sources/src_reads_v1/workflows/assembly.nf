include { READS_ASSEMBLE; READS_ANNOTATE; ASSEMBLY_COLLECT } from '../modules/local/reads'
workflow ASSEMBLY {
    take: clean_units; databases; code; code_digest
    main:
    READS_ASSEMBLE(clean_units,code,code_digest)
    branches = READS_ASSEMBLE.out.units.branch {
        annotate: params.enable_phispy && it[1]=='bacterial_isolate'
        plain: true
    }
    READS_ANNOTATE(branches.annotate,databases,code,code_digest)
    combined = branches.plain.mix(READS_ANNOTATE.out.units).map { id,mat,dir -> dir }
    ASSEMBLY_COLLECT(combined.collect().ifEmpty([]),code,code_digest)
    emit: assemblies = ASSEMBLY_COLLECT.out.assemblies
}
