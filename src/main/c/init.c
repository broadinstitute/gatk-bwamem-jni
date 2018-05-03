#include "init.h"

jfieldID peStatClass_failedID;
jfieldID peStatClass_lowID;
jfieldID peStatClass_highID;
jfieldID peStatClass_averageID;
jfieldID peStatClass_stdID;


#define FAIL_IF_NULL(x) if (!(x)) return JNI_ERR;

jint JNI_OnLoad(JavaVM* vm, void* reserved) {

    JNIEnv* env;
    jclass peStatClass;
    if ((*vm)->GetEnv(vm, &env, JNI_VERSION_1_8) != JNI_OK) {
         return JNI_ERR;
    }

    FAIL_IF_NULL(peStatClass = (*env)->FindClass(env, "org/broadinstitute/hellbender/utils/bwa/BwaMemPairEndStats"));
    FAIL_IF_NULL(peStatClass_failedID = (*env)->GetFieldID(env, peStatClass, "failed", "Z"));
    FAIL_IF_NULL(peStatClass_lowID = (*env)->GetFieldID(env, peStatClass, "low", "I"));
    FAIL_IF_NULL(peStatClass_highID = (*env)->GetFieldID(env, peStatClass, "high", "I"));
    FAIL_IF_NULL(peStatClass_averageID = (*env)->GetFieldID(env, peStatClass, "average", "D"));
    FAIL_IF_NULL(peStatClass_stdID = (*env)->GetFieldID(env, peStatClass, "std", "D"));
    (*env)->DeleteLocalRef(env, peStatClass);

    return JNI_VERSION_1_8;
}

