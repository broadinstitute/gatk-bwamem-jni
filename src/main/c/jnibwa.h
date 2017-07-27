/*
 * jnibwa.h
 */

#ifndef JNIBWA_H_
#define JNIBWA_H_

#include "bwa/bwamem.h"


int jnibwa_createReferenceIndex( char const* refFileName, char const* indexPrefix, char const* algoName);
int jnibwa_createIndexFile( char const* refName, char const* imgSuffix );
bwaidx_t* jnibwa_openIndex( int fd );
int jnibwa_destroyIndex( bwaidx_t* pIdx );
void* jnibwa_getRefContigNames( bwaidx_t* pIdx, size_t* pBufSize );
void* jnibwa_createAlignments( bwaidx_t* pIdx, mem_opt_t* pOpts, char* pSeq, size_t* pBufSize);

#endif /* JNIBWA_H_ */
