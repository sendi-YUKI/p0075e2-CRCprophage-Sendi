process ASSOCIATION_MODELS {
    tag "${spec.simpleName}"
    cpus 2
    memory '4 GB'
    maxForks 1
    container params.association_container
    publishDir "${params.outdir}", mode:'copy', overwrite:true
    input:
    path long_table
    path metadata
    path spec
    path schema
    path python_helper
    path r_helper
    output:
    path 'association', emit: result
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    python ${python_helper} --long ${long_table} --metadata ${metadata} --spec ${spec} \
      --schema ${schema} --r-script ${r_helper} --outdir association
    """
}
