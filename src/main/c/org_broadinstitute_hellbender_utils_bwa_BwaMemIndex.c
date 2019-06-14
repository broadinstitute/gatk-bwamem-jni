#include <jni.h>
#include <fcntl.h>
#include <stdlib.h>
#include "jnibwa.h"
#include "bwa/bwa_commit.h"

/*
 * Implementation of native routines declared in BwaMemIndex.java.
 * These fucntions take care of interacting with the JNI environment, primarily marshalling and unmarshalling arguments.
 * It delegates to functions in jnibwa.c to do the interface with the BWA-MEM API.
 */

static char* jstring_to_chars( JNIEnv* env, jstring in ) {
    const char* tmp = (*env)->GetStringUTFChars(env, in, 0);
    char* res = strdup(tmp);
    (*env)->ReleaseStringUTFChars(env, in, tmp);
    return res;
}

static void throwErrorMessage( JNIEnv* env, char* message ) {
    jclass exceptionClass = (*env)->FindClass(env, "org/broadinstitute/hellbender/utils/bwa/BwaMemException");
    if ( exceptionClass ) {
        (*env)->ThrowNew(env, exceptionClass, message);
        free(message);
    }
}


JNIEXPORT void JNICALL Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_indexReference(
                        JNIEnv* env, jclass cls, jstring jReferenceFileName, jstring jIndexPrefix, jint jAlgo ) {
    char* refName = jstring_to_chars(env, jReferenceFileName);
    char* indexPrefix = jstring_to_chars(env, jIndexPrefix);
    char* errMsg = 0;
    jnibwa_indexReference(refName, indexPrefix, jAlgo, &errMsg);
    free(refName);
    free(indexPrefix);
    if ( errMsg ) throwErrorMessage(env, errMsg);
}

JNIEXPORT void JNICALL Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_createIndexImageFile(
                        JNIEnv* env, jclass cls, jstring referencePrefix, jstring imageFileName ) {
    char* refName = jstring_to_chars(env, referencePrefix);
    char* imgName = jstring_to_chars(env, imageFileName);
    char* errMsg = 0;
    jnibwa_createIndexFile(refName, imgName, &errMsg);
    free(refName);
    free(imgName);
    if ( errMsg ) throwErrorMessage(env, errMsg);
}

JNIEXPORT jlong JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_openIndex( JNIEnv* env, jclass cls, jstring memImgFilename ) {
    char* imageName = jstring_to_chars(env, memImgFilename);
    char* errMsg = 0;
    jlong address = (jlong)jnibwa_openIndex(imageName, &errMsg);
    free(imageName);
    if ( errMsg ) throwErrorMessage(env, errMsg);
    return address;
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_destroyIndex( JNIEnv* env, jclass cls, jlong idxAddr ) {
    if ( !idxAddr ) {
        throwErrorMessage(env, strdup("null index address"));
    } else {
        char* errMsg = 0;
        jnibwa_destroyIndex((bwaidx_t*)idxAddr, &errMsg);
        if ( errMsg ) throwErrorMessage(env, errMsg);
    }
}

JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_createDefaultOptions( JNIEnv* env, jclass cls ) {
    size_t defaultsSize = sizeof(mem_opt_t)+sizeof(mem_pestat_t);
    void* defaults = malloc(defaultsSize);
    void* opts = mem_opt_init();
    memcpy(defaults, opts, sizeof(mem_opt_t));
    free(opts);
    jobject result = (*env)->NewDirectByteBuffer(env, defaults, defaultsSize);
    if ( !result ) throwErrorMessage(env, strdup("unable to create ByteBuffer for default options"));
    return result;
}

// returns a ByteBuffer with the reference contig names
// returned ByteBuffer has:
//   a 32-bit int giving the number of contig names
//   for each contig name,
//     a 32-bit int giving the length of the name
//     that many bytes of name
JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_getRefContigNames( JNIEnv* env, jclass cls, jlong idxAddr ) {
    jobject namesBuf = 0;
    if ( !idxAddr ) {
        throwErrorMessage(env, strdup("null index address"));
    } else {
        size_t bufSize = 0;
        char* errMsg = 0;
        void* bufMem = jnibwa_getRefContigNames((bwaidx_t*)idxAddr, &bufSize, &errMsg);
        if ( errMsg ) {
            free(bufMem);
            throwErrorMessage(env, errMsg);
        } else {
            namesBuf = (*env)->NewDirectByteBuffer(env, bufMem, bufSize);
            if ( !namesBuf ) {
                free(bufMem);
                throwErrorMessage(env, strdup("unable to create ByteBuffer for reference contig names"));
            }
        }
    }
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
    if ( !idxAddr ) {
        throwErrorMessage(env, strdup("null index address"));
        return 0;
    }
    mem_opt_t* pOpts = (*env)->GetDirectBufferAddress(env, optsBuf);
    if ( !pOpts ) {
        throwErrorMessage(env, strdup("can't get address for opts ByteBuffer"));
        return 0;
    }
    char* pSeq = (*env)->GetDirectBufferAddress(env, seqsBuf);
    if ( !pSeq ) {
        throwErrorMessage(env, strdup("can't get address for seqs ByteBuffer"));
        return 0;
    }
    size_t bufSize = 0;
    void* bufMem;
    mem_pestat_t* inniesStats = (mem_pestat_t*)(pOpts+1);
    char* errMsg = 0;
    if ( !((char*)inniesStats)[-1] ) {
        bufMem = jnibwa_createAlignments((bwaidx_t*)idxAddr, pOpts, 0, pSeq, &bufSize, &errMsg);
    } else {
        mem_pestat_t peStats[4];
        memset(peStats, 4*sizeof(mem_pestat_t), 0);
        peStats[0].failed = peStats[2].failed = peStats[3].failed = 1;
        memcpy(peStats+1, inniesStats, sizeof(mem_pestat_t));
        bufMem = jnibwa_createAlignments((bwaidx_t*)idxAddr, pOpts, peStats, pSeq, &bufSize, &errMsg);
    }
    if ( errMsg ) {
        free(bufMem);
        throwErrorMessage(env, errMsg);
        return 0;
    }
    jobject alnBuf = (*env)->NewDirectByteBuffer(env, bufMem, bufSize);
    if ( !alnBuf ) {
        free(bufMem);
        throwErrorMessage(env, strdup("can't create ByteBuffer for alignments"));
    }
    return alnBuf;
}

JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_createByteBuffer( JNIEnv* env, jclass cls, jint bufSize ) {
    void* bufMem = malloc(bufSize);
    if ( !bufMem ) {
        throwErrorMessage(env, strdup("can't allocate memory for ByteBuffer"));
        return 0;
    }
    jobject result = (*env)->NewDirectByteBuffer(env, bufMem, bufSize);
    if ( !result ) throwErrorMessage(env, strdup("can't create ByteBuffer"));
    return result;
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_destroyByteBuffer( JNIEnv* env, jclass cls, jobject alnBuf ) {
    void* buf = (*env)->GetDirectBufferAddress(env, alnBuf);
    if ( !buf ) throwErrorMessage(env, strdup("can't get ByteBuffer address"));
    else free(buf);
}

JNIEXPORT jstring JNICALL
Java_org_broadinstitute_hellbender_utils_bwa_BwaMemIndex_getVersion( JNIEnv* env, jclass cls ) {
    jstring result = (*env)->NewStringUTF(env, BWA_COMMIT);
    if ( !result ) throwErrorMessage(env, strdup("can't create version string"));
    return result;
}
