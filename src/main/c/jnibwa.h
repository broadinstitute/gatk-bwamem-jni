/*
 * jnibwa.h
 */

#ifndef JNIBWA_H_
#define JNIBWA_H_

#include "bwa/bwamem.h"


void jnibwa_indexReference( char const* refFileName, char const* indexPrefix, int algo, char** pErrMsg );
void jnibwa_createIndexFile( char const* refName, char const* imgSuffix, char** pErrMsg );
bwaidx_t* jnibwa_openIndex( char const* imgName, char** pErrMsg );
void jnibwa_destroyIndex( bwaidx_t* pIdx, char** pErrMsg );
void* jnibwa_getRefContigNames( bwaidx_t* pIdx, size_t* pBufSize, char** pErrMsg );
void* jnibwa_createAlignments( bwaidx_t* pIdx, mem_opt_t* pOpts, mem_pestat_t* peStats, char* pSeq, size_t* pBufSize, char** pErrMsg );

#endif /* JNIBWA_H_ */
