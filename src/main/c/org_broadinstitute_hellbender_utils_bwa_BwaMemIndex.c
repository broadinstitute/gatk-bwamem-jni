#include <jni.h>
#include <fcntl.h>
#include <stdlib.h>
#include "jnibwa.h"
#include "bwa/bwa_commit.h"


char * jstring_to_const_char(JNIEnv* env, jstring in) {
    const char* tmp = (*env)->GetStringUTFChars(env, in, 0);
    char* res = strdup(tmp);
    (*env)->ReleaseStringUTFChars(env, in, tmp);
    return res;
}


JNIEXPORT jboolean JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_createReferenceIndex( JNIEnv* env, jclass cls, jstring jReferenceFileName, jstring jIndexPrefix, jstring jAlgoName) {

    extern int bwt_idx_build(const char*, const char*, const char*);

	char *referenceFileName = jstring_to_const_char(env, jReferenceFileName);
	char *indexPrefix = jstring_to_const_char(env, jIndexPrefix);
	char *algoName = jstring_to_const_char(env, jAlgoName);
	jboolean res = !bwt_idx_build( referenceFileName, indexPrefix, algoName);
	free(referenceFileName); free(indexPrefix); free(algoName);
	return res;
}

JNIEXPORT jboolean JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_createIndexImageFile( JNIEnv* env, jclass cls, jstring referencePrefix, jstring imageFileName ) {
	char *refName = jstring_to_const_char(env, referencePrefix);
	char *imgName = jstring_to_const_char(env, imageFileName);
	jboolean res = !jnibwa_createIndexFile( refName, imgName );
	free(refName); free(imgName);
	return res;
}

JNIEXPORT jlong JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_openIndex( JNIEnv* env, jclass cls, jstring memImgFilename ) {
	char *fname = jstring_to_const_char(env, memImgFilename);
	int fd = open(fname, O_RDONLY);
	free(fname);
	if ( fd == -1 ) return 0;
	return (jlong)jnibwa_openIndex(fd);
}

JNIEXPORT jint JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_destroyIndex( JNIEnv* env, jclass cls, jlong idxAddr ) {
	if ( !idxAddr ) return 0;
	return jnibwa_destroyIndex((bwaidx_t*)idxAddr);
}

JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_createDefaultOptions( JNIEnv* env, jclass cls ) {
	return (*env)->NewDirectByteBuffer(env, mem_opt_init(), sizeof(mem_opt_t));
}

// returns a ByteBuffer with the reference contig names
// returned ByteBuffer has:
//   a 32-bit int giving the number of contig names
//   for each contig name,
//     a 32-bit int giving the length of the name
//     that many bytes of name
JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_getRefContigNames( JNIEnv* env, jclass cls, jlong idxAddr ) {
	if ( !idxAddr ) return 0;
	size_t bufSize = 0;
	void* bufMem = jnibwa_getRefContigNames((bwaidx_t*)idxAddr, &bufSize);
	jobject namesBuf = (*env)->NewDirectByteBuffer(env, bufMem, bufSize);
	if ( !namesBuf ) free(bufMem);
	return namesBuf;
}

// we accept a ByteBuffer that contains:
//   a 32-bit integer count of the number of sequences to follow
//   each sequence is just a regular old C string (8-bit characters, null terminated) giving the bases in the sequence
// the idxAddr is what you got from the createIndex method above
// the optsBuf argument is a mem_opt_t structure wrapped by a ByteBuffer (from createDefaultOptions method)
// we return a ByteBuffer that contains:
// for each sequence,
//   a 32-bit integer count of the number of alignments that follow
//   for each alignment, a flattened BAM-like pseudo-structure like this:
/*
typedef struct {
	int32_t flag_mapQ; // flag<<16 | mapQ (the flag value is a SAM-formatted flag)

	// these fields are present only if the read is mapped, i.e., flag&4==0
	int32_t refID; // reference id
	int32_t pos; // reference starting position (0-based)
	int32_t NM; // value for NM tag (number of mismatches)
	int32_t AS; // value for AS tag (alignment score)
	int32_t XS; // value for bwa-specific XS tag (suboptimal alignment score)
	int32_t nCigar; // nCigarOps
	int32_t cigarOp[nCigarOps]; // len<<4 | op
	int32_t nMDchars; // length of MD tag
	char    MDtag[(nMDchars+3)&~3]; // value for MD tag (space allocation rounded up to stay on int32_t boundary)
	int32_t nXAchars; // length of XA tag
	char    XAtag[(nXAchars+3)&~3]; // value for XA tag (space rounded up to int32_t boundary)

	// these fields are present only if the read is paired and the mate is mapped, i.e., flag&9==1
	int32_t mateRefID; // mate's reference id
	int32_t matePos;   // mate's reference starting position (0-based)
	int32_t tlen;      // inferred template length
} Alignment;
*/
JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_createAlignments(
				JNIEnv* env, jclass cls, jobject seqsBuf, jlong idxAddr, jobject optsBuf ) {
	bwaidx_t* pIdx = (bwaidx_t*)idxAddr;
	mem_opt_t* pOpts = (*env)->GetDirectBufferAddress(env, optsBuf);
	char* pSeq = (*env)->GetDirectBufferAddress(env, seqsBuf);
	size_t bufSize = 0;
	void* bufMem = jnibwa_createAlignments(pIdx, pOpts, pSeq, &bufSize);
	jobject alnBuf = (*env)->NewDirectByteBuffer(env, bufMem, bufSize);
	if ( !alnBuf ) free(bufMem);
	return alnBuf;
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_destroyByteBuffer( JNIEnv* env, jclass cls, jobject alnBuf ) {
	free((*env)->GetDirectBufferAddress(env, alnBuf));
}

JNIEXPORT jstring JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_getVersion( JNIEnv* env, jclass cls ) {
	return (*env)->NewStringUTF(env, BWA_COMMIT);
}
