#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "incl/klib/khash.h"
#include "incl/klib/kseq.h"
#include "incl/klib/kvec.h"
#include "vcf.h"
#include <zlib.h>

/*
 * sim_genome
 *
 * Jeremy Wang
 * 20191106
 *
 * As fast as possible, simulate a genome/genotype given a reference
 * sequence and VCF with alleles and frequencies
*/

// have to reorder params to make this work with kseq
int fileread(FILE* f, char* buffer, int size) {
  return fread(buffer, 1, size, f);
}

// init kseq struct
KSEQ_INIT(gzFile, gzread)

// creates string:[array of uint8] hash
// to map chrom names to sequences
KHASH_MAP_INIT_STR(refSeq, char*);

// creates string:uint32 hash
// to map chrom names to lengths
KHASH_MAP_INIT_STR(refLen, uint32_t);

// creates string:int hash
// to map chrom names to IDs (tid)
KHASH_MAP_INIT_STR(refId, int);

int main(int argc, char *argv[]) {
  srand(time(0));

  if(argc < 2) {
    fprintf(stderr, "Usage: sg <reference FASTA> <VCF>\n");
    fprintf(stderr, "Not enough arguments.\n");
    return 1;
  }
  char *ref_fasta = argv[1];
  char *vcf_file = argv[2];

  khint_t bin, bin2, bin3; // hash bin (result of kh_put)
  int absent;

  // load ref FASTA file
  //
  khash_t(refSeq) *ref = kh_init(refSeq);
  khash_t(refLen) *rlen = kh_init(refLen);
  khash_t(refId) *rid = kh_init(refId);

  gzFile gzfp;
  kseq_t *seq;
  int l;

  gzfp = gzopen(ref_fasta, "r");
  if(!gzfp) {
    fprintf(stderr, "File '%s' not found\n", ref_fasta);
    return 1;
  }
  fprintf(stderr, "Reading fasta file: %s\n", ref_fasta);
  seq = kseq_init(gzfp);

  int refid = 0;
  char* dup;
  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //printf("Reading %s (%i bp).\n", seq->name.s, l);

    dup = malloc(sizeof(char) * (strlen(seq->name.s) + 1));
    dup[strlen(seq->name.s)] = '\0';
    memcpy(dup, seq->name.s, sizeof(char) * strlen(seq->name.s));

    // seq array
    bin = kh_put(refSeq, ref, dup, &absent);
    // copy the seq read from kseq to a new heap here - this is pretty fast and the easiest way to implement right now (see kseq.h)
    kh_val(ref, bin) = malloc(sizeof(char)*l);
    memcpy(kh_val(ref, bin), seq->seq.s, sizeof(char)*l);

    // sequence length
    bin = kh_put(refLen, rlen, dup, &absent);
    kh_val(rlen, bin) = l;

    // ref ID (indexed order in FASTA)
    bin = kh_put(refId, rid, dup, &absent);
    kh_val(rid, bin) = refid++;
  }

  gzclose(gzfp);
  kseq_destroy(seq);

  int i, j;

  // load VCF file
  //
  fprintf(stderr, "Reading VCF file '%s'\n", vcf_file);
  vcf_file_t vcf = vcf_init(vcf_file);
  if(vcf.fp == NULL) {
    fprintf(stderr, "ERROR: Failed reading VCF file '%s', check that it exists\n", vcf_file);
    return 1;
  }

  // check that Chroms match what's in our ref FASTA
  char modchrom[100];
  uint64_t found_bp = 0;
  uint64_t unfound_bp = 0;
  for(i = 0; i < vcf.header.num_sequences; i++) {
    fprintf(stderr, "%s\n", vcf.header.sequences[i]);
    bin = kh_get(refSeq, ref, vcf.header.sequences[i]);
    if(bin == kh_end(ref)) { // absent
      sprintf(modchrom, "chr%s", vcf.header.sequences[i]);
      bin = kh_get(refSeq, ref, modchrom);
      if(bin == kh_end(ref)) { // absent
        fprintf(stderr, "WARNING: chromosome '%s' in VCF not found in ref FASTA\n", vcf.header.sequences[i]);
        unfound_bp = unfound_bp + vcf.header.sequence_lengths[i];
        continue;
      } else {
        // put the names WITHOUT chr into the hash, reusing the VCF header names and pointing to the same seq
        bin2 = kh_put(refSeq, ref, vcf.header.sequences[i], &absent);
        kh_val(ref, bin2) = kh_val(ref, bin);
        bin3 = kh_get(refLen, rlen, modchrom);
        fprintf(stderr, "%s 10m-10m100: ", vcf.header.sequences[i]);
        for(j = 10000000; j < 10000100; j++) {
          fprintf(stderr, "%c", kh_val(ref, bin2)[j]);
        }
        fprintf(stderr, "\n");
      }
    } else {
      bin3 = kh_get(refLen, rlen, vcf.header.sequences[i]);
    }
    if(vcf.header.sequence_lengths[i] != kh_val(rlen, bin3)) {
      fprintf(stderr, "ERROR: chromosome '%s' sizes don't match: %u in VCF, %u in FASTA\n", vcf.header.sequence_lengths[i], kh_val(rlen, bin3));
      unfound_bp = unfound_bp + vcf.header.sequence_lengths[i];
      return 1;
    }
    found_bp = found_bp + vcf.header.sequence_lengths[i];
  }
  if((double)found_bp / (found_bp + unfound_bp) < 0.9) {
    fprintf(stderr, "ERROR: we seem to have matched < 90\% of the expected bp from the VCF with the reference FASTA, maybe you have the wrong reference\n");
    return 1;
  }

  fprintf(stderr, "Processing variants...\n");
  uint32_t n_snps = 0;
  uint32_t variants_added = 0;
  vcf_line_t *entry = vcf_read_line(&vcf);
  while(entry != NULL) {
    n_snps++;
    if(strlen(entry->ref) == 1 && strlen(entry->alt) == 1 && ((double)rand()) / RAND_MAX < entry->af) {
      bin = kh_get(refSeq, ref, entry->chrom);
      if(bin != kh_end(ref)) {
        if(kh_val(ref, bin)[entry->pos-1] != entry->ref[0] && kh_val(ref, bin)[entry->pos-1]-32 != entry->ref[0]) { // may be upper or lower case in ref
          fprintf(stderr, "ERROR: %s -> %s at %s:%d (but we found a %c)\n", entry->ref, entry->alt, entry->chrom, entry->pos, kh_val(ref, bin)[entry->pos-1]);
        }
        //fprintf(stderr, "Adding %s -> %s at %s:%d\n", entry->ref, entry->alt, entry->chrom, entry->pos);
        variants_added++;
      }
    }
    vcf_destroy_line(entry);

    entry = vcf_read_line(&vcf);
    if(n_snps % 1000000 == 0) {
      fprintf(stderr, "\r["); // beginning of ghetto progress bar
      for(i = 0; i < n_snps/1000000; i++)
        fprintf(stderr, "#");
      for(; i < 84; i++)
        fprintf(stderr, " ");
      fprintf(stderr, "]");
    }
  }
  fprintf(stderr, "\n");
  vcf_close(&vcf);

  fprintf(stderr, "%u SNPs evaluated\n", n_snps);
  fprintf(stderr, "%u variants added\n", variants_added);

  // free the hashes, but be aware half the keys in the ref hash are already free because we reused them from the VCF

  return 0;
}

