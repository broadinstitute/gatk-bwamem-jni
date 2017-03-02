#!/usr/bin/env bash

#Upload a snapshot build of the current HEAD with both .dylib and .so
#This can only be run from within Broad

#set -v
set -e

SERVER=gsa6.broadinstitute.org

LIB_PATH="src/main/c"

echo "building local files"
./gradlew clean build

echo "building remote files"
vectorLib=$( ssh $SERVER 'bash -s' < scripts/build_native_lib_in_clean_repo.sh $( git rev-parse HEAD) gatk-bwamem-jni libbwa $LIB_PATH )

echo "result is at $vectorLib"

echo "copying from ${SERVER}:${vectorLib} to $LIB_PATH"
scp ${SERVER}:${vectorLib} $LIB_PATH

exit 0
