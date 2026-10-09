nextflow.enable.dsl=2
include { CRC_PHAGE_CORE } from './workflows/crc_phage_core'

workflow {
    if (!params.a1_result) error 'Required --a1_result: successful A0+A1 output directory'
    a1 = file(params.a1_result, checkIfExists:true)
    master = a1.resolve('prophage_master.tsv')
    if (!master.exists()) error 'A1 prophage_master.tsv missing'
    statuses = a1.resolve('sample_status.tsv').readLines().findAll { it.trim() }
    headers = statuses[0].split('\t', -1).toList()
    idx = headers.indexOf('genome_id')
    if (idx < 0 || statuses.size() < 2) error 'A1 sample_status.tsv missing genome IDs'
    bankFile = file(params.genbanks, checkIfExists:true)
    bankRows = bankFile.readLines().findAll { it.trim() }
    bankHeader = bankRows[0].split('\t', -1).toList()
    if (!bankHeader.contains('genome_id') || !bankHeader.contains('genbank')) error 'genbanks TSV needs genome_id and genbank'
    banks = [:]
    bankSeen = [] as Set
    bankRows.drop(1).each { line ->
        cells = line.split('\t', -1)
        gid = cells[bankHeader.indexOf('genome_id')]
        value = cells[bankHeader.indexOf('genbank')]
        if (bankSeen.contains(gid)) error "Duplicate GenBank ID: ${gid}"
        bankSeen.add(gid)
        if (value) {
            bank = java.nio.file.Paths.get(value)
            if (!bank.isAbsolute()) bank = bankFile.parent.resolve(value)
            bank = file(bank)
            if (!bank.exists()) error "GenBank file missing: ${bank}"
            banks[gid] = bank
        }
    }
    samples = statuses.drop(1).collect { line ->
        id = line.split('\t', -1)[idx]
        prepared = a1.resolve("samples/${id}/input/prepared")
        if (!prepared.resolve('genome.fna').exists()) error "Missing A1 prepared FASTA: ${id}"
        tuple([id:id], prepared)
    }.sort { left, right -> left[0].id <=> right[0].id }
    if (!samples.collect { it[0].id }.containsAll(bankSeen)) error 'GenBank manifest contains IDs not in A1 samples'
    if (samples.collect { it[0].id }.unique().size() != samples.size()) error 'Duplicate A1 sample IDs'
    checkm2Db = file(params.checkm2_db, checkIfExists:true)
    functionsDb = file(params.functions_db, checkIfExists:true)
    refs = file(params.reference_dir, checkIfExists:true)
    metadata = file(params.metadata, checkIfExists:true)
    CRC_PHAGE_CORE(Channel.fromList(samples), a1, banks, checkm2Db, functionsDb, refs, metadata)
}

workflow.onComplete {
    def statusDir = new File("${params.outdir}/pipeline_info")
    statusDir.mkdirs()
    new File(statusDir,'workflow_status.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        success:workflow.success, exit_status:workflow.exitStatus, session_id:workflow.sessionId.toString(),
        run_name:workflow.runName, current_target:'A2_A3_A4', completed:workflow.complete.toString()
    ])) + '\n'
}
