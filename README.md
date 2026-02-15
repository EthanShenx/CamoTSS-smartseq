# CamoTSS: Analysis of Alternative Transcription Start Sites

CamoTSS is a tool for analyzing alternative transcription start sites (TSS) for cellular phenotypes and regulatory patterns from 5' scRNA-seq data.

## Note

Hi there, my github account did not notify me when there are issue.
So if you are in a hurry, you can email me. ruiyan@connect.hku.hk.
I check email every day.

## Installation

You can install from this GitHub repository for latest (often development)
version by following command line

```bash
pip install -U git+https://github.com/StatBiomed/CamoTSS
```

In either case, add `--user` if you don't have the write permission for your
Python environment.

## Quick start

### Download test file

You can download test file from [figshare](https://figshare.com/projects/CamoTSS/184603).

Here, you can download some large file include genome.fa, possorted_genome_bam_filtered.bam.

## Run CamoTSS

Here are three modes in CamoTSS: **TC+CTSS**, **TC** and **CTSS**.

When you run **TC+CTSS** mode, you will get TC result and then get the CTSS result based on the TC.

When you run **TC** mode, you will only get the TSS cluster result.

The **TC+CTSS** and **TC** mode have the same required files.

The --outdir is the only required parameter for **CTSS** mode. But the outdir should include output of TC.

If you want to run **CTSS** mode, you must based on the output of TC.

CamoTSS supports two sequencing platforms: **10x** (default) and **smartseq5**.

### 10x Genomics Data Example:

For the remaining modes, you can check the [full documentation](https://camotss.readthedocs.io/en/latest/run_CamoTSS.html).

```bash
#!/bin/bash
gtfFile=$download/Homo_sapiens.GRCh38.105.chr_test.gtf
fastaFile=$download/genome.fa
bamFile=$download/possorted_genome_bam_filtered.bam
cellbarcodeFile=$download/cellbarcode_to_CamoTSS

CamoTSS --gtf $gtfFile --refFasta $fastaFile --bam $bamFile -c $cellbarcodeFile -o CamoTSS_out --mode TC+CTSS
```

### Smart-seq5 Data Example:

For Smart-seq5 data, each cell is represented by a separate BAM file. You can run CamoTSS with the --platform smartseq5 option:

```bash
#!/bin/bash
gtfFile=$download/Homo_sapiens.GRCh38.105.chr_test.gtf
fastaFile=$download/genome.fa
bamDir=/path/to/smartseq5/bams/  # Directory containing individual BAM files (one per cell)
outDir=./CamoTSS_smartseq5_out

CamoTSS --platform smartseq5 --gtf $gtfFile --refFasta $fastaFile --bam_dir $bamDir -o $outDir --mode TC+CTSS
```

Alternatively, you can provide a list of BAM files:

```bash
#!/bin/bash
gtfFile=$download/Homo_sapiens.GRCh38.105.chr_test.gtf
fastaFile=$download/genome.fa
bamList=/path/to/bam_list.txt  # File containing paths to BAM files (one per line)
outDir=./CamoTSS_smartseq5_out

CamoTSS --platform smartseq5 --gtf $gtfFile --refFasta $fastaFile --bam_list $bamList -o $outDir --mode TC+CTSS
```

### Smart-seq5 Options:

* `--platform {10x,smartseq5}` - Specify sequencing platform (default: 10x)
* `--bam_list <file>` - File containing list of BAM files for smartseq5 mode (one BAM per line)
* `--bam_dir <dir>` - Directory containing BAM files for smartseq5 mode (each BAM represents one cell)
* `--cell_id_from {filename,tsv}` - How to determine cell ID for smartseq5: from BAM filename (default) or from TSV mapping
* `--cell_map <cells.tsv>` - TSV file mapping sample names to cell IDs for smartseq5 mode
* `--dedup {umi,coord,fragment,none}` - Deduplication method: umi (for 10x), coord/fragment (for smartseq5), none. Default depends on platform.
* `--min_mapq <int>` - Minimum mapping quality for reads (default: 20)
* `--tss_read {read1,read2}` - Which mate contains the 5' transcript sequence used for TSS calling. Default: read1 for 10x, read2 for smartseq5.

**Note**
You should use the same reference gtf file and reference fasta file as that you used during alignment. In other words, if you run alignment by using cellranger, then the gtf file and fasta file should located in the refdata-gex-GRCh38-2020-A/fasta/genome.fa and refdata-gex-GRCh38-2020-A/genes/genes.gtf.

## Alternative TSS or CTSS detecting

To identify alternative TSS usage or alternative CTSS usage, Brie2 (Huang & Sanguinetti, 2021) is recommend to be used.

Here, we provide an example exploiting BRIE2 to detect alterntive TSS/CTSS usage.

You can check it in our [manual](https://camotss.readthedocs.io/en/latest/runBRIE.html).

## Detailed Manual

The full manual is [here](https://camotss.readthedocs.io/en/latest/index.html), including:

[Preprocess](https://camotss.readthedocs.io/en/latest/preprocess.html)

[Run CamoTSS](https://camotss.readthedocs.io/en/latest/run_CamoTSS.html)

[Detect alternative TSS/CTSS](https://camotss.readthedocs.io/en/latest/runBRIE.html)

## Reference

Hou, R., Hon, CC. & Huang, Y. CamoTSS: analysis of alternative transcription start sites for cellular phenotypes and regulatory patterns from 5' scRNA-seq data. Nat Commun 14, 7240 (2023). https://doi.org/10.1038/s41467-023-42636-1
