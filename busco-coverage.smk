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

def mosdepth_regions():
    reg = config["regions"]
    if reg is None:
        return ""
    else:
        return f"--by {reg}"

#  ------------------------------------------------------------- TARGET
rule all:
    input: config["prefix"] + ".mosdepth.global.dist.txt"



# --------------------------------------------------------------------
# RULES
# --------------------------------------------------------------------


localrules: cleanup_reference
rule cleanup_reference:
    input: config["reference"]
    output: temp(config["prefix"] + ".ref.fa")
    shell: "seqkit replace -p '.+' -r 'NODE_{{nr}}' {input} > {output}"


#  --------------------------------------------------- busco diamond db
localrules: diamond_db
rule diamond_db:
    input: config["busco"]
    output: temp(config["prefix"] + ".dmnd")
    params: config["prefix"]
    shell: "diamond makedb --db {params} --in {input}"


#  ---------------------------------------------------- index reference
rule index:
    input: rules.cleanup_reference.output
    output: temp(multiext(rules.cleanup_reference.output[0], ".fai", ".sa", ".amb", ".pac", ".bwt", ".ann"))
    log: ".log/" + config["prefix"] + ".idx.log"
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
        cram=temp(config["prefix"] + ".cram"),
    log: ".log/" + config["prefix"] + ".map.log"
    resources:
        runtime = runtimer(60),
        mem_mb = mem_gb(8)
    threads: 8
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.R1} {input.R2} 2>> {log} \
        | samtools sort -@ {threads} --reference {input.ref} - 2>> {log} \
        | samtools view -@ {threads} -C -T {input.ref} - > {output}
        """


#  --------------------------------------------------------- index cram
localrules: index_cram
rule index_cram:
    input:
        cram=rules.map_reads.output.cram,
        ref=rules.cleanup_reference.output,
        refi=rules.index.output
    output:
        crami=temp(config["prefix"] + ".cram.crai")
    threads: 4
    shell: "samtools index {input.cram}"


# -------------------------------------------- find busco in scaffolds
localrules: find_buscos
rule find_buscos:
    input:
        ref=rules.cleanup_reference.output,
        db=rules.diamond_db.output
    output:
        blast=config["prefix"] + ".busco.blast"
    params:
        prefix=config["prefix"],
        minlength=config["blast_minlength"]
    resources:
        runtime=runtimer(2),
        mem_mb=mem_gb(4)
    shell:
        "diamond blastx --db {params.prefix} --query {input.ref} 2> /dev/null "
        " > {output} "

#  ---------------------------------------------- convert blastx to BED
localrules: blast_to_bed
rule blast_to_bed:
    input: rules.find_buscos.output.blast
    output: temp(config["prefix"] + ".busco.hsp.bed")
    shell:
        "cut -f1,9,10 {input} "
        "| bedtools sort -i - "
        "| awk '{{if($3-$2 > {params.minlength}) print}}' "
        " > {output}"

#  --------------------------------------------------- compute coverage
rule mosdepth:
    input:
        cram=rules.map_reads.output.cram,
        crami=rules.index_cram.output.crami,
        ref=rules.cleanup_reference.output,
        refi=rules.index.output,
        regions=rules.blast_to_bed.output.bed
    output:
        dist=config["prefix"] + ".mosdepth.global.dist.txt"
    params:
        prefix=config["prefix"],
    log: ".log/" + config["prefix"] + ".mosdepth.log"
    resources:
        mem_mb = mem_gb(8),
        runtime = runtimer(20)
    threads: 4
    shell:
        "mosdepth -x -n "
        " --by {input.regions} "
        " --threads {threads} "
        " --fasta {input.ref} "
        " {params.prefix} {input.cram} 2> {log}"
