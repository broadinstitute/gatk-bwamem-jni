#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
//#include <time.h>
#include <unistd.h>
#include <zlib.h>
#include <errno.h>

#include "bwa/kstring.h"
#include "bwa/bntseq.h"
#include "bwa/bwt.h"
#include "bwa/utils.h"
#include "bwa/rle.h"
#include "bwa/rope.h"

int bwt_idx_build(const char *fasta, const char *prefix, const char *algo_type_str);

int bwt_idx_build(const char *fasta, const char *prefix, const char *algo_type_str)
{
	extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);
	extern bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is);

	char *str, *str2, *str3;
	int algo_type = 0;
// Commented out clocking and debug print outs, but kept them just in
//	case they are needed at some point.
//	clock_t t;
	int64_t l_pac;

	if (strcmp(algo_type_str, "rb2") == 0) algo_type = 1;
	else if (strcmp(algo_type_str, "is") == 0) algo_type = 3;
	else if (strcmp(algo_type_str, "auto") == 0) algo_type = 0;
	else err_fatal(__func__, "unknown algorithm: '%s'.", algo_type_str);

	str  = (char*)calloc(strlen(prefix) + 10, 1);
	str2 = (char*)calloc(strlen(prefix) + 10, 1);
	str3 = (char*)calloc(strlen(prefix) + 10, 1);

	{ // nucleotide indexing
		gzFile fp = xzopen(fasta, "r");
//		t = clock();
//		fprintf(stderr, "[bwa_index] Pack FASTA... ");
		l_pac = bns_fasta2bntseq(fp, prefix, 0);
//		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		err_gzclose(fp);
	}
	if (algo_type == 0) algo_type = l_pac > 50000000? 1 : 3; // set the algorithm for generating BWT
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".pac");
		strcpy(str2, prefix); strcat(str2, ".bwt");
//		t = clock();
//		fprintf(stderr, "[bwa_index] Construct BWT for the packed sequence...\n");
		bwt = bwt_pac2bwt(str, algo_type == 3);
		bwt_dump_bwt(str2, bwt);
		bwt_destroy(bwt);
//		fprintf(stderr, "[bwa_index] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".bwt");
//		t = clock();
//		fprintf(stderr, "[bwa_index] Update BWT... ");
		bwt = bwt_restore_bwt(str);
		bwt_bwtupdate_core(bwt);
		bwt_dump_bwt(str, bwt);
		bwt_destroy(bwt);
//		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	{
		gzFile fp = xzopen(fasta, "r");
//		t = clock();
//		fprintf(stderr, "[bwa_index] Pack forward-only FASTA... ");
		l_pac = bns_fasta2bntseq(fp, prefix, 1);
//		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		err_gzclose(fp);
	}
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".bwt");
		strcpy(str3, prefix); strcat(str3, ".sa");
//		t = clock();
//		fprintf(stderr, "[bwa_index] Construct SA from BWT and Occ... ");
		bwt = bwt_restore_bwt(str);
		bwt_cal_sa(bwt, 32);
		bwt_dump_sa(str3, bwt);
		bwt_destroy(bwt);
//		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	free(str3); free(str2); free(str);
	return 0;
}
