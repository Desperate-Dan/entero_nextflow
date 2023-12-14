process ampliClean {
  
  publishDir path: "${params.out_dir}/${barcode}/ampli_clean", mode: 'copy'

  input:
    tuple val(barcode), path(binned_reads)
    path refs
    path bed
    val min
    val max
    
    
  output:
    path("${barcode}.*.fastq.gz"), emit: reads, optional: true
    path "log.txt"

  script:
    """
    ampli_clean -f ${binned_reads} -r ${refs} -o ${barcode} -b ${bed} --min ${min} --max ${max} --map-only -s --fastq --log
    """
}

process medaka_hap_var {

  publishDir path: "${params.out_dir}/${barcode}/medaka", mode: 'copy'

  input:
    path (input_reads)
    path (refs)
    val (medaka_model)

  output:
    path "${barcode}.${vir}.consensus.fasta"

  script:
    barcode = input_reads.name.toString().tokenize('.').get(0)
    vir = input_reads.name.toString().tokenize('.').get(1)
    """
    grep -A 1 ${vir} ${refs} > ${vir}.ref.fasta
    medaka_haploid_variant -i ${input_reads} -r ${vir}.ref.fasta -m ${medaka_model} -f -x
    medaka stitch medaka/consensus_probs.hdf ${vir}.ref.fasta ${barcode}.${vir}.consensus.fasta
    """
}
//These lines for fastq dir parsing are taken from rmcolq's workflow https://github.com/rmcolq/pantheon
EXTENSIONS = ["fastq", "fastq.gz", "fq", "fq.gz"]

ArrayList get_fq_files_in_dir(Path dir) {
    return EXTENSIONS.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
}

workflow {
//Define input channels  
  ref_ch = file("${params.refs}")
  bed_ch = file("${params.bed}")
  schemes_dir_ch = file("${params.schemes_dir}")
  min_ch = Channel.value("${params.min}")
  max_ch = Channel.value("${params.max}")
  med_mod_ch = Channel.value("${params.medaka_model}")
//These lines for fastq dir parsing are taken from rmcolq's workflow https://github.com/rmcolq/pantheon
  run_dir = file("${params.fastq}", type: "dir", checkIfExists:true)
  barcode_input = Channel.fromPath("${run_dir}/*", type: "dir", checkIfExists:true, maxDepth:1).map { [it.baseName, get_fq_files_in_dir(it)]}

//Run the processes
  ampliClean(barcode_input, ref_ch, bed_ch, min_ch, max_ch)
  medaka_hap_var(ampliClean.out.reads.flatten(), ref_ch, med_mod_ch)
}