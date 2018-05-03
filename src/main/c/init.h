#include <jni.h>

#ifndef INIT_H
#define INIT_H

extern jfieldID peStatClass_failedID;
extern jfieldID peStatClass_lowID;
extern jfieldID peStatClass_highID;
extern jfieldID peStatClass_averageID;
extern jfieldID peStatClass_stdID;


jint JNI_OnLoad(JavaVM* vm, void* reserved);

#endif