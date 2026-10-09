nextflow.enable.dsl=2
include { ASSEMBLY } from './workflows/assembly'
workflow {
    if (!params.clean_reads_manifest || !params.outdir) error 'Require --clean_reads_manifest and --outdir'
    def manifest=file(params.clean_reads_manifest,checkIfExists:true)
    def data=new groovy.json.JsonSlurper().parse(manifest.toFile())
    if(data.schema_version!=1 || !(data.units instanceof List)) error 'Invalid or empty clean manifest'
    def seenUnits = new HashSet()
    def resolveRead = { value ->
        def raw = value.toString()
        file(new File(raw).isAbsolute() ? raw : "${manifest.parent}/${raw}", checkIfExists:true).toAbsolutePath().normalize()
    }
    def units=data.units.collect { original ->
        def u = new LinkedHashMap(original)
        if (!(u.unit_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]{0,100}/) || !seenUnits.add(u.unit_id)) error 'Invalid/duplicate assembly unit_id'
        if (!(u.material_type in ['bacterial_isolate','bulk_metagenome','vlp'])) error 'Unsupported assembly material_type'
        def originals = ([u.reads_1,u.reads_2]+(u.singletons ?: [])).findAll{it}
        if (!u.reads_1 || !(u.checksums instanceof Map) || originals.toSet()!=u.checksums.keySet() || originals.toSet().size()!=originals.size()) error 'Clean unit read paths/checksums disagree'
        def resolved = originals.collectEntries { value -> [(value):resolveRead(value).toString()] }
        u.reads_1=resolved[u.reads_1]
        u.reads_2=u.reads_2 ? resolved[u.reads_2] : null
        u.singletons=(u.singletons ?: []).collect { resolved[it] }
        u.checksums=u.checksums.collectEntries { path, digest -> [(resolved[path]):digest] }
        def unitB64=groovy.json.JsonOutput.toJson(u).bytes.encodeBase64().toString()
        tuple(u.unit_id,u.material_type,unitB64,originals.collect { file(resolved[it],checkIfExists:true) })
    }
    def codeNames = ['reads_io.py','reads_process.py','reads_assembly.py','reads_publish.py']
    def codeContent = codeNames.collect { n -> n + ':' + java.security.MessageDigest.getInstance('SHA-256').digest(file("${projectDir}/bin/${n}",checkIfExists:true).bytes).encodeHex().toString() }.join('\n')
    def codeDigest = java.security.MessageDigest.getInstance('SHA-256').digest(codeContent.bytes).encodeHex().toString()
    ASSEMBLY(Channel.fromList(units),file(params.database_manifest,checkIfExists:true),file("${projectDir}/bin"),codeDigest)
}
