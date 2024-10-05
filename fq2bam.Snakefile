from glob import glob
configfile: "conf.yaml"
ref=config["ref"]
trimdir=config["trimdir"]
samps=["2367","2416","2462","AGU_3339_12","Brisas-1","CAB_71_PL3","CAB_76_PL3","CAB_77_PL5","CATIE_1000","CCAT1119_EET544_A1","CCAT1858_EET547","CCAT4675_EET575","CCAT4688_EET576_A1","CCAT4998_EET577","CRU_89","CUR3_G35_A1","CUR3_G37_A6","CUR3_G38_A8","EB2237_A1","EET103_Borde","EET446_A1","EET451_A1","EET462_A1","EET58_A1","EET95_A2","EET_103","EET_395","EET_400","EET_59","K82","KA2","LCT46_A1","LCT_EEN_46","LP_3_40","LP_4_48","LX_32","PBC123","PMF_20"]
rule all:
  input:
    expand(os.path.join(f"{trimdir}","{sample}.bam"),sample=samps),
  threads: 4
rule trim:
  input:
    fw="{sample}_1.fastq.gz",
    rv="{sample}_2.fastq.gz",
    idx=config["ref"]
  output:
    trim1=expand(os.path.join(f"{trimdir}","{sample}_1_val_1.fq.gz"),sample=samps),
    trim2=expand(os.path.join(f"{trimdir}","{sample}_2_val_2.fq.gz"),sample=samps)
  threads: 4
  shell:
    """
    trim_galore --output_dir trim --quality 28 --illumina --phred33 --stringency 6 --length 60 --paired {input.fw} {input.rv}
    """
rule aln:
  input:
    ref=config["ref"],
    fw=expand(os.path.join(f"{trimdir}","{sample}_1_val_1.fq.gz"),sample=samps),
    rv=expand(os.path.join(f"{trimdir}","{sample}_2_val_2.fq.gz"),sample=samps)
  output:
    fw=expand(os.path.join(f"{trimdir}","{sample}_1.sai"),sample=samps),
    rv=expand(os.path.join(f"{trimdir}","{sample}_2.sai"),sample=samps),
    sam=expand(os.path.join(f"{trimdir}","{sample}.sam"),sample=samps)
  params:
    readgroup=r"\@RG\tID:{sample}\tSM:{sample}"
  shell:
    """
    bwa aln -n 0.06 -M 6 {input.ref} {input.fw} > {output.fw}
    bwa aln -n 0.06 -M 6 {input.ref} {input.rv} > {output.rv}
    bwa sampe -r '{params.readgroup}' -t 4 {input.ref} {output.fw} {output.rv} {input.fw} {input.rv} > {output.sam}
    """
rule viewbam:
  input:
    sam="{sample}.sam",
    ref=config["ref"]
  output:
    bam="{sample}.bam"
  shell:
    """
    samtools view -ubhSt {input.ref} {input.sam} -o {output.bam}
    """
