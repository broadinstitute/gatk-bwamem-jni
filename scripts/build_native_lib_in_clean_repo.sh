#!/usr/bin/env bash

#Helper script to be called from build_both_dylib_and_dylib.sh
#
#checkout and clone a new copy of gatk in a tmpdir
#build the .so file and return its location
#this can only be run on broad machines
#
# usage build_native_lib_in_clean_repo.sh <commit hash> <git-project> <library name> <path to lib>

set -v

COMMIT=$1
PROJECT=$2
LIB_NAME=$3
LIB_PATH=$4


source /broad/tools/scripts/useuse
reuse Git-2.11

set -e

export TMPDIR="/broad/hptmp"

echoerr() { echo "$@" 1>&2; }

PROJECT_TMP_DIR=`mktemp --tmpdir -d 2>/dev/null || mktemp -d -t 'mytmpdir'`

echoerr "Moving to $PROJECT_TMP_DIR"
cd "$PROJECT_TMP_DIR"

echoerr "cloning broadinstitute/${PROJECT}"
GIT_LFS_SKIP_SMUDGE=1 git clone git@github.com:broadinstitute/${PROJECT}.git 1>2

echoerr "Moving to $PROJECT"
cd $PROJECT

echoerr "Checking out ${commit}"
GIT_LFS_SKIP_SMUDGE=1 git checkout -f "$COMMIT" 1>2

echoerr "performing gradle build"
./gradlew build 1>2 --no-daemon

LIB=$(pwd)/$(find $LIB_PATH -name "${LIB_NAME}*" ! -name "*.a" )
echoerr "created $LIB"

echo $LIB
exit 0
