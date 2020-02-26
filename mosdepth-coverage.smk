configfile: "config.yaml"


def runtimer(time):
    def sk(w, attempt):
        return int(attempt**attempt * time)
    return sk

def add_cpu(ncpu):
    def cpu(w, attempt):
        return int(attempt * ncpu)
    return cpu

def mem_gb(gb):
    mb = int(gb * 1024)
    def f(w, attempt):
        return int(attempt * mb)
    return f


#  ----------------------------------------------- cleanup fasta header
localrules: cleanup_reference
rule cleanup_reference:
    input: config["reference"]
    output: config["prefix"] + ".ref.fa"
    shell: "seqkit replace -p '.+' -r 'NODE_{{nr}}' {input} > {output}"


#  ---------------------------------------------------- index reference
rule index:
    input: rules.cleanup_reference.output
    output: temp(multiext(rules.cleanup_reference.output[0], ".fai", ".sa", ".amb", ".pac", ".bwt"))
    log: config["prefix"] + ".idx.log"
    resources:
        runtime = runtimer(10),
        mem_mb = mem_gb(4)
    threads: 1
    shell:
        """
        bwa index {input} 2>> {log}
        samtools faidx {input} 2>> {log}
        """

#  ---------------------------------------------------------- map reads
rule map_reads:
    input:
        ref=rules.cleanup_reference.output,
        refi=rules.index.output,
        R1=config["R1"],
        R2=config["R2"]
    output:
        bam=temp(config["prefix"] + ".bam"),
    log: config["prefix"] + ".map.log"
    resources:
        runtime = runtimer(60),
        mem_mb = mem_gb(8)
    threads: 8
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.R1} {input.R2} 2>> {log} \
        | samtools sort -@ {threads} -O BAM --reference {input.ref} -o {output.bam} - 2>> {log}
        """

#  ------------------------------------------------ convert bam to cram
localrules: bam_to_cram
rule bam_to_cram:
    input:
        bam=rules.map_reads.output.bam,
        ref=rules.cleanup_reference.output
    output:
        cram=config["prefix"] + ".cram"
    threads: 4
    shell:
        "samtools view -C -T {input.ref} {input.bam} > {output}"

#  --------------------------------------------------------- index cram
localrules: index_cram
rule index_cram:
    input:
        cram=rules.bam_to_cram.output.cram,
        ref=rules.cleanup_reference.output,
        refi=rules.index.output
    output:
        crami=temp(config["prefix"] + ".cram.crai")
    threads: 4
    shell: "samtools index {input.cram}"

#  -------------------------------------------------- compute coverage
rule mosdepth:
    input:
        cram=rules.bam_to_cram.output.cram,
        ref=rules.cleanup_reference.output,
        refi=rules.index.output,
        crami=rules.index_cram.output.crami
    output:
        dist=config["prefix"] + ".mosdepth.global.dist.txt"
    params:
        prefix=config["prefix"],
        regions=mosdepth_regions()
    log: config["prefix"] + ".mosdepth.log"
    resources:
        mem_mb = mem_gb(8),
        runtime = runtimer(20)
    threads: 4
    shell:
        "mosdepth --fast-mode --no-per-base "
        " {params.regions} "
        " --threads {threads} --fasta {input.ref}"
        " {params.prefix} {input.cram} 2> {log}"
