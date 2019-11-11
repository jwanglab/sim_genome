Simulate a genome with a random genotype derived from a VCF file with allele frequencies.

Right now, it only adds single-nucleotide variants.


Installation
------------

    git clone http://github.com/jwanglab/sim_genome
    cd sim_genome
    mkdir incl
    cd incl
    git clone http://github.com/attractivechaos/klib
    cd ..
    make
    ./sg -h


Usage
-----

    ./sg ref.fasta snps.vcf > output.fasta


Example
-------

This example should run in a couple of minutes (after data is downloaded) and use just over 3GB of RAM.

    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    ./sg hg19.fa.gz ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz > hg19_1kgenomes_genotype.fasta

By default, the random seed used to generate each random genotype is the current epoch and is printed to stderr.
