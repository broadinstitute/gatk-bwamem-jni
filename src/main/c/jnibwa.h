/*
 * jnibwa.h
 */

#ifndef JNIBWA_H_
#define JNIBWA_H_

#include "bwa/bwamem.h"

char* jnibwa_createIndexFile( char const* refName, char const* imgSuffix );
char* jnibwa_openIndex( char const* imgName, int ignoreVersion, int compareCRC );
char* jnibwa_destroyIndex( bwaidx_t* pIdx );
void* jnibwa_getRefContigNames( bwaidx_t* pIdx, size_t* pBufSize );
void* jnibwa_createAlignments( bwaidx_t* pIdx, mem_opt_t* pOpts, char* pSeq, size_t* pBufSize);

#endif /* JNIBWA_H_ */
