/*
 * jnibwa.c
 */

#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <errno.h>
#include "jnibwa.h"
#include "bwa/kstring.h"

static inline void kput32( int32_t val, kstring_t* str ) {
	kputsn((char*)&val, sizeof(int32_t), str);
}

static inline int cigarRefLen( int nCigar, uint32_t const* pCigar )
{
	uint32_t const* pEnd = pCigar + nCigar;
	int len = 0;
	while ( pCigar != pEnd ) {
		uint32_t lenOp = *pCigar++;
		int op = lenOp & 0xf;
		if ( !op || op == 2 )
			len += lenOp >> 4;
	}
	return len;
}

static void fmt_BAMish(mem_opt_t const* opt, bntseq_t const* bns, kstring_t *str, bseq1_t *s, int n, mem_aln_t const* list, int which, mem_aln_t const* p, mem_aln_t const* m) {
	if ( !which ) {
		size_t nInts = 12; // for space planning, assume mapped, unpaired reads with 3 cigar ops and an 8 character MD
		if ( p->flag & 0x1 ) nInts += 3; // if paired, add in enough space for mate info
		ks_resize(str, (n*nInts+1)*sizeof(int32_t));
		kput32(n, str);
	}
	int32_t flag_mapQ = p->flag;
	if ( p->flag & 0x10000 ) flag_mapQ |= 0x100;
	flag_mapQ = (flag_mapQ << 16) | (p->mapq & 0xff);
	kput32(flag_mapQ, str);
	if ( !(p->flag & 0x4) ) {
		kput32(p->rid, str);
		kput32(p->pos, str);
		kput32(p->NM, str);
		kput32(p->score, str);
		kput32(p->sub, str);
		int32_t nCig = p->n_cigar;
		kput32(nCig, str);
		uint32_t* pCig = p->cigar;
		while ( nCig-- ) {
			uint32_t lenOp = *pCig++;
			// op is encoded as MIDSH in a mem_aln_t, but as MIDNSH in a BAM
			if ( (lenOp & 0xf) > 2 ) ++lenOp;
			kput32(lenOp, str);
		}
		int32_t nMD = p->n_cigar ? strlen((char*)pCig) : 0;
		kput32(nMD, str);
		if ( nMD )
			kputsn((char*)pCig, (nMD+3)&~3, str);
		int32_t nXA = p->XA ? strlen(p->XA) : 0;
		kput32(nXA, str);
		if ( nXA ) {
			kputsn(p->XA, (nXA+3)&~3, str);
		}
	}
	if ( (p->flag & 0x9) == 1 ) {
		kput32(m->rid, str);
		kput32(m->pos, str);
		if ( (p->flag & 0x4) || p->rid != m->rid ) kput32(0, str);
		// the next two lines represent my interpretation of the SAM spec
		// else if ( p->pos < m->pos ) kput32(m->pos+cigarRefLen(m->n_cigar, m->cigar)-p->pos, str);
		// else kput32(m->pos-p->pos-cigarRefLen(p->n_cigar, p->cigar), str);
		// but BWA does something else which is very odd in the case of outies,
		// but is faithfully reproduced below
		else {
			// get 5' end of + reads, and 3' end of - reads
			long p0 = p->pos;
			if ( p->is_rev ) p0 += cigarRefLen(p->n_cigar, p->cigar) - 1;
			long m0 = m->pos;
			if ( m->is_rev ) m0 += cigarRefLen(m->n_cigar, m->cigar) - 1;
			kput32(m0 - p0 + (p0 > m0 ? -1 : p0 < m0 ? 1 : 0), str);
		}
	}
}

static size_t bufLen( int32_t* pBuf ) {
	size_t totLen = 1;
	size_t nAligns = *pBuf++;
	while ( nAligns-- ) {
		int32_t flag = *pBuf++ >> 16;
		totLen += 1; // for flags
		if ( !(flag & 0x4) ) {
			totLen += 8; // refId, pos, NM, AS, XS, nCigOps, nMDchars, nXAchars
			pBuf += 5;
			int32_t nCig = *pBuf++;
			totLen += nCig;
			pBuf += nCig;
			int32_t nMD = (*pBuf++ + 3) >> 2;
			totLen += nMD;
			pBuf += nMD;
			int32_t nXA = (*pBuf++ + 3) >> 2;
			totLen += nXA;
			pBuf += nXA;
		}
		if ( (flag & 0x9) == 1 ) {
			totLen += 3; // mate rid, mate pos, tlen
			pBuf += 3;
		}
	}
	return totLen;
}

int jnibwa_createIndexFile( char const* refName, char const* imgName ) {
  	bwaidx_t* pIdx = bwa_idx_load(refName, BWA_IDX_ALL);
  	bwa_idx2mem(pIdx);
  	int fd = open(imgName, O_WRONLY|O_CREAT|O_TRUNC, 0644);
  	if ( fd == -1 ) {
  		printf("Failed to open %s for writing: %s\n", imgName, strerror(errno));
  		return 2;
  	}
  	size_t len = pIdx->l_mem;
  	uint8_t* buf = pIdx->mem;
  	while ( len ) {
  		size_t toWrite = len;
  		if ( toWrite > (1L<<30) ) toWrite = 1L<<30;
  		if ( write(fd,buf,toWrite) != toWrite ) {
  	  		printf("Failed to write %s: %s\n", imgName, strerror(errno));
  	  		return 2;
  		}
  		buf += toWrite;
  		len -= toWrite;
  	}
  	if ( close(fd) != 0 ) {
  		printf("Failed to close %s: %s\n", imgName, strerror(errno));
  		return 2;
  	}
	bwa_idx_destroy(pIdx);
	return 0;
}

bwaidx_t* jnibwa_openIndex( int fd ) {
	struct stat statBuf;
	if ( fstat(fd, &statBuf) == -1 ) return 0;
	uint8_t* mem = mmap(0, statBuf.st_size, PROT_READ, MAP_SHARED, fd, 0);
	close(fd);
	bwaidx_t* pIdx = calloc(1, sizeof(bwaidx_t));
	bwa_mem2idx(statBuf.st_size, mem, pIdx);
	pIdx->is_shm = 1;
	mem_fmt_fnc = &fmt_BAMish;
	bwa_verbose = 0;
	return pIdx;
}

int jnibwa_destroyIndex( bwaidx_t* pIdx ) {
	void* pMem = pIdx->mem;
	size_t memLen = pIdx->l_mem;
	bwa_idx_destroy(pIdx);
	return munmap(pMem, memLen);
}

void* jnibwa_getRefContigNames( bwaidx_t* pIdx, size_t* pBufSize ) {
	int nRefContigs = pIdx->bns->n_seqs;
	bntann1_t* pAnnoBeg = pIdx->bns->anns;
	bntann1_t* pAnnoEnd = pAnnoBeg + nRefContigs;
	bntann1_t* pAnno;
	int bufSize = 4 + 4*nRefContigs; // for the ints that describe the number of contigs, and the length of each name
	for ( pAnno = pAnnoBeg; pAnno != pAnnoEnd; ++pAnno ) {
		bufSize += strlen(pAnno->name)+1;
	}
	char* bufMem = malloc(bufSize);
	*(int32_t*)bufMem = nRefContigs;
	char* pBuf = bufMem + sizeof(int32_t);
	for ( pAnno = pAnnoBeg; pAnno != pAnnoEnd; ++pAnno ) {
		size_t len = strlen(pAnno->name);
		*(int32_t*)pBuf = len;
		pBuf += sizeof(int32_t);
		memcpy(pBuf, pAnno->name, len);
		pBuf += len;
	}
	*pBufSize = bufSize;
	return bufMem;
}

void* jnibwa_createAlignments( bwaidx_t* pIdx, mem_opt_t* pOpts, char* pSeq, size_t* pBufSize) {
	char c = 0;
	char* emptyString = &c;
	uint32_t nSeqs = *(uint32_t*)pSeq;
	pSeq += sizeof(uint32_t);
	bseq1_t* pSeq1Beg = calloc(nSeqs, sizeof(bseq1_t));
	bseq1_t* pSeq1End = pSeq1Beg+nSeqs;
	bseq1_t* pSeq1;
	for ( pSeq1 = pSeq1Beg; pSeq1 != pSeq1End; ++pSeq1 ) {
		size_t seqLen = strlen(pSeq);
		pSeq1->l_seq = seqLen;
		pSeq1->seq = pSeq;
		pSeq1->name = emptyString;
		pSeq1->id = pSeq1-pSeq1Beg;
		pSeq += seqLen + 1;
	}

	mem_process_seqs(pOpts, pIdx->bwt, pIdx->bns, pIdx->pac, 0, nSeqs, pSeq1Beg, 0);

	size_t nInts = 0;
	for ( pSeq1 = pSeq1Beg; pSeq1 != pSeq1End; ++pSeq1 ) {
		if ( pSeq1->sam ) nInts += bufLen((int32_t*)pSeq1->sam);
	}
	int32_t* resultsBeg = malloc(nInts*sizeof(int32_t));
	int32_t* pResults = resultsBeg;
	for ( pSeq1 = pSeq1Beg; pSeq1 != pSeq1End; ++pSeq1 ) {
		int32_t* pBuf = (int32_t*)pSeq1->sam;
		if ( pBuf ) {
			size_t len = bufLen(pBuf);
			memcpy(pResults, pBuf, len*sizeof(int32_t));
			free(pBuf);
			pResults += len;
		}
	}
	free(pSeq1Beg);

	*pBufSize = nInts*sizeof(int32_t);
	return resultsBeg;
}
