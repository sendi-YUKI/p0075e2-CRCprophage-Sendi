process CORE_QC_COLLECT {
    cpus 1
    memory '1 GB'
    container 'python:3.12.13-bookworm@sha256:3cd9086bdb30f7c9bc08a3fa621d9842e0d3f6f9291aeb4677e0547817c10b12'
    stageInMode 'copy'
    publishDir "${params.outdir}/host_qc", mode:'copy'
    input:
    path qc_dirs, stageAs:'qc??/*'
    path core_helper, stageAs:'core_genomes.py'
    path phage_helper, stageAs:'phageflow.py'
    output:
    path 'host_qc', emit: qc
    script:
    def inputs = (qc_dirs instanceof List ? qc_dirs : [qc_dirs]).collect { "'${it}/genome_qc.tsv'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    mkdir -p host_qc
    python3 core_genomes.py deduplicate --inputs ${inputs} --outdir host_qc
    """
}

process CORE_BLAST {
    cpus 6
    memory '4 GB'
    container 'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759@sha256:7f0f18604de90429273e086528b050ad4db0e8bb2eff5a63576afc38e222d27b'
    publishDir "${params.outdir}/raw/votu", mode:'copy'
    input:
    path fasta, stageAs:'candidates.fna'
    output:
    path 'blast', emit: blast
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    mkdir -p blast
    blastn -version > blast/version.txt
    if [ -s '${fasta}' ]; then
      nseq=\$(grep -c '^>' '${fasta}')
      makeblastdb -in '${fasta}' -dbtype nucl -out blast/database > blast/makeblastdb.log 2>&1
      blastn -task blastn -query '${fasta}' -db blast/database -evalue 1e-5 \\
        -max_target_seqs "\$nseq" -num_threads ${task.cpus} \\
        -outfmt '6 std qlen slen qseq sseq' -out blast/alignments.tsv \\
        > blast/blastn.stdout.log 2> blast/blastn.stderr.log
      echo completed > blast/status.txt
    else
      touch blast/alignments.tsv
      echo zero_candidates > blast/status.txt
    fi
    """
}

process CORE_VOTU {
    cpus 1
    memory '2 GB'
    container 'python:3.12.13-bookworm@sha256:3cd9086bdb30f7c9bc08a3fa621d9842e0d3f6f9291aeb4677e0547817c10b12'
    publishDir "${params.outdir}/catalog", mode:'copy'
    input:
    path master, stageAs:'baseline_master.tsv'
    path fasta, stageAs:'candidates.fna'
    path blast
    path qc
    path element_class, stageAs:'element_class.tsv'
    path helper, stageAs:'core_votu.py'
    val policy
    output:
    path 'catalog', emit: catalog
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 core_votu.py --master '${master}' --fasta '${fasta}' \\
      --blast '${blast}/alignments.tsv' --qc '${qc}/genome_qc.tsv' \\
      --ani-percent ${policy.ani_percent} --af-shorter-percent ${policy.af_shorter_percent} --element-class '${element_class}' --require-element-evidence --outdir catalog
    """
}
