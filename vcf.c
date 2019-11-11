#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "vcf.h"

vcf_header_t parse_header(FILE* vcf_fp) {
  vcf_header_t header;
  // 1024 maximum ref sequences
  header.num_sequences = 0;
  header.sequences = malloc(sizeof(char*) * 1024);
  header.sequence_lengths = malloc(sizeof(size_t) * 1024);

  char line[10000]; // maximum line size is 1024 chars
  const char* eq = "=";
  const char* comma = ",";
  const char* semic = ";";
  char *key, *val;
  
  // get lines as long as they start with '#' (header)
  int c;
  int i, j;
  uint64_t l;
  while(1) {
    c = fgetc(vcf_fp); // get first character
    ungetc(c, vcf_fp); // QUICK! put it back! (basically just a peek)
    if((char)c == '#') { // header line
      if(fgets(line, sizeof line, vcf_fp) != NULL) {
        char* ln = malloc(sizeof(char) * strlen(line));
        memcpy(ln, line, sizeof(char) * strlen(line));
        ln[strlen(line)-1] = '\0'; // these lines will always end in a newline, so just truncate it
        if(ln[1] != '#') { // this is the column header - we'll assume we know what it says
          continue;
        }
        for(i = 0; i < strlen(line); i++) {
          if(ln[i] == '=') {
            ln[i] = '\0';
            if(strcmp(ln+2, "contig") == 0) {
              // ##contig=<ID=1,assembly=b37,length=249250621>
              //fprintf(stderr, "found contig: %s\n", ln+i+1);
              key = ln+i+2; // starts right after the '<'
              for(i = i + 2; i < strlen(line); i++) {
                if(ln[i] == '=') {
                  ln[i] = '\0';
                  val = ln+i+1;
                } else if(ln[i] == ',' || ln[i] == '>') {
                  ln[i] = '\0';
                  if(strcmp(key, "ID") == 0) {
                    header.sequences[header.num_sequences] = malloc((strlen(val)+1) * sizeof(char));
                    memcpy(header.sequences[header.num_sequences], val, strlen(val)+1);
                  }
                  if(strcmp(key, "length") == 0) {
                    sscanf(val, "%llu", &l); // covert val (str) to uint64_t
                    header.sequence_lengths[header.num_sequences] = l;
                  }
                  key = ln+i+1;
                }
              }
              header.num_sequences++;
            }
            if(strcmp(ln+2, "fileformat") == 0) {
              // ##fileformat=VCFv4.1
              //fprintf(stderr, "VCF format: %s\n", ln+i+1);
              header.version = malloc((strlen(ln+i+1)+1) * sizeof(char));
              memcpy(header.version, ln+i+1, sizeof(char) * (strlen(ln+i+1)+1));
            }
            break;
          }
        }
        free(ln);
      }
    } else { // non-header line
      break;
    }
  }

  return header;
}

vcf_line_t *vcf_read_line(vcf_file_t* vcf) {
  //const char* vcf_format = "%s\t%d\t%d\t%s\t%d\t%c";

  vcf_line_t *line = malloc(sizeof(vcf_line_t));
  
  vcf->cur_row++;

  // parse vcf line by tokenizing - allows arbitrary fields after vcf-6
  char* s = NULL; // this MUST be initialized to NULL or getline with complain
  const char delim[3] = "\t\n";
  size_t slen;
  int32_t nread = getline(&s, &slen, vcf->fp);
  if(nread <= 0) return NULL;
  char* token;
  char* key;
  char* val;
  token = strtok(s, delim);
  int i, j;
  // #CHROM  POS ID  REF ALT QUAL  FILTER  INFO
  for(i = 0; token != NULL; i++) {
    //fprintf(stderr, "token %d: %s\n", i, token);
    switch(i) {
      case 0:
        line->chrom = malloc((strlen(token)+1) * sizeof(char));
        strcpy(line->chrom, token);
        line->chrom[strlen(token)] = '\0'; // manually add null-terminator
        break;
      case 1:
        line->pos = atoi(token);
        break;
      case 2:
        break;
      case 3:
        line->ref = malloc((strlen(token)+1) * sizeof(char));
        strcpy(line->ref, token);
        line->ref[strlen(token)] = '\0'; // manually add null-terminator
        break;
      case 4:
        line->alt = malloc((strlen(token)+1) * sizeof(char));
        strcpy(line->alt, token);
        line->alt[strlen(token)] = '\0'; // manually add null-terminator
        break;
      case 5:
      case 6:
        break;
      case 7: // INFO
        // semicolon-delimited KEY=value
        // AC=2130;AF=0.425319;AN=5008;NS=2504;DP=103152;EAS_AF=0.3363;AMR_AF=0.3602;AFR_AF=0.4909;EUR_AF=0.4056;SAS_AF=0.4949;AA=|||unknown(NO_COVERAGE);VT=INDEL
        key = token;
        int l = strlen(token);
        for(j = 0; j < l+1; j++) {
          if(token[j] == '=') {
            token[j] = '\0';
            val = token+j+1;
          } else if(token[j] == ';' || token[j] == '\0') {
            token[j] = '\0';
            if(strcmp(key, "AF") == 0) {
              sscanf(val, "%f", &(line->af));
              //fprintf(stderr, "AF: %f\n", line->af);
            }
            key = token+j+1;
          }
        }
        break;
    }
    token = strtok(NULL, delim);
  }
  free(s);

  return line;
}

int vcf_destroy_line(vcf_line_t* l) {
  free(l->chrom);
  free(l->ref);
  free(l->alt);
  free(l);
}

vcf_file_t vcf_init(char* f) {

  // open vcf file from path
  FILE *vcf_fp = fopen(f, "r");
  if( vcf_fp == NULL ) {
    fprintf(stderr, "Error reading vcf file '%s'\n", f);
  }

  vcf_header_t header;
  fprintf(stderr, "Parsing header\n");
  header = parse_header(vcf_fp);

  fprintf(stderr, "VCF version %s\n", header.version);
  fprintf(stderr, "%d sequences found\n", header.num_sequences);

  vcf_file_t vcf = {vcf_fp, header, 0};

  return vcf;
}

int vcf_close(vcf_file_t* vcf) {
  fclose(vcf->fp);
  if(vcf->header.sequences != NULL) {
    free(vcf->header.sequences);
  }
  if(vcf->header.sequence_lengths != NULL) {
    free(vcf->header.sequence_lengths);
  }
}
