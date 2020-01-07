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

// creates string:int hash
// to map chrom names to indices (into refs and rlens)
KHASH_MAP_INIT_STR(refName, uint32_t);

int main(int argc, char *argv[]) {

  if(argc < 2) {
    fprintf(stderr, "Usage: sg <reference FASTA> <VCF>\n");
    fprintf(stderr, "Not enough arguments.\n");
    return 1;
  }
  char *ref_fasta = argv[1];
  char *vcf_file = argv[2];
  int seed = argc > 2 ? atoi(argv[3]) : 0;
  if(seed == 0) {
    time_t now = time(0);
    fprintf(stderr, "Random seed (time): %ld\n", now);
    srand(now);
  } else {
    fprintf(stderr, "Random seed: %ld\n", seed);
    srand(seed);
  }

  khint_t bin, bin2, bin3; // hash bin (result of kh_put)
  int absent;

  // load ref FASTA file
  //
  khash_t(refName) *ref_lookup = kh_init(refName);

  kvec_t(char*) refs;
  kv_init(refs);

  kvec_t(uint32_t) rlens;
  kv_init(rlens);

  kvec_t(char*) rnames;
  kv_init(rnames);

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

  uint32_t refids = 0;
  uint32_t rid;
  char* dup;
  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //fprintf(stderr, "Reading %s (%i bp).\n", seq->name.s, l);

    dup = malloc(sizeof(char) * (strlen(seq->name.s) + 1));
    dup[strlen(seq->name.s)] = '\0';
    memcpy(dup, seq->name.s, sizeof(char) * strlen(seq->name.s));

    // seq name lookup
    bin = kh_put(refName, ref_lookup, dup, &absent);
    kh_val(ref_lookup, bin) = refids++;
    kv_push(char*, rnames, dup);

    // copy the seq read from kseq to a new heap here - this is pretty fast and the easiest way to implement right now (see kseq.h)
    char* s = malloc(sizeof(char)*(l+1));
    memcpy(s, seq->seq.s, sizeof(char)*l);
    s[l] = '\0';
    kv_push(char*, refs, s);

    // sequence length
    kv_push(uint32_t, rlens, l);
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
    bin = kh_get(refName, ref_lookup, vcf.header.sequences[i]);
    if(bin == kh_end(ref_lookup)) { // absent
      sprintf(modchrom, "chr%s", vcf.header.sequences[i]);
      bin = kh_get(refName, ref_lookup, modchrom);
      if(bin == kh_end(ref_lookup)) { // absent
        fprintf(stderr, "WARNING: chromosome '%s' in VCF not found in ref FASTA\n", vcf.header.sequences[i]);
        unfound_bp = unfound_bp + vcf.header.sequence_lengths[i];
        continue;
      } else {
        rid = kh_val(ref_lookup, bin);
        // put the names WITHOUT chr into the hash, reusing the VCF header names
        bin2 = kh_put(refName, ref_lookup, vcf.header.sequences[i], &absent);
        kh_val(ref_lookup, bin2) = rid;
      }
    } else {
      rid = kh_val(ref_lookup, bin);
    }
    if(vcf.header.sequence_lengths[i] != kv_A(rlens, rid)) {
      fprintf(stderr, "ERROR: chromosome '%s' sizes don't match: %u in VCF, %u in FASTA\n", vcf.header.sequence_lengths[i], kv_A(rlens, rid));
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
  uint32_t prev_pos = 0;
  char prev_ref;
  vcf_line_t *entry = vcf_read_line(&vcf);
  while(entry != NULL) {
    n_snps++;
    if(strlen(entry->ref) == 1 && strlen(entry->alt) == 1 && ((double)rand()) / RAND_MAX < entry->af) {
      bin = kh_get(refName, ref_lookup, entry->chrom);
      if(bin != kh_end(ref_lookup)) {
        rid = kh_val(ref_lookup, bin);
        // check if our ref seq agrees with the VCF
        if(kv_A(refs, rid)[entry->pos-1] != entry->ref[0] && kv_A(refs, rid)[entry->pos-1]-32 != entry->ref[0]) { // may be upper or lower case in ref
          // note: sometimes there are two variants at the same site and it may have been changed already...
          if(prev_pos != entry->pos || prev_ref != entry->ref[0]) {
            fprintf(stderr, "ERROR: %s -> %s at %s (rid %u):%d (but we found %c)\n", entry->ref, entry->alt, entry->chrom, rid, entry->pos, kv_A(refs, rid)[entry->pos-1]);
            break;
          }
        }
        //fprintf(stderr, "Adding %s -> %s at %s:%d\n", entry->ref, entry->alt, entry->chrom, entry->pos);
        kv_A(refs, rid)[entry->pos-1] = entry->alt[0];
        variants_added++;
      }
      prev_pos = entry->pos;
      prev_ref = entry->ref[0];
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

  // print new seqs
  fprintf(stderr, "Writing new FASTA file (to stdout)...\n");
  for(rid = 0; rid < refids; rid++) {
    fprintf(stdout, ">%s\n", kv_A(rnames, rid));
    for(i = 0; i < kv_A(rlens, rid); i++) {
      fprintf(stdout, "%c", kv_A(refs, rid)[i]);
      if(i%80 == 79) {
        fprintf(stdout, "\n");
      }
    }
    if(i%80 != 80) {
      fprintf(stdout, "\n");
    }
  }

  // free the hashes, but be aware half the keys in the ref hash are already free because we reused them from the VCF

  return 0;
}

