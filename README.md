[![Maven Central](https://maven-badges.herokuapp.com/maven-central/org.broadinstitute/gatk-bwamem-jni/badge.svg)](https://maven-badges.herokuapp.com/maven-central/org.broadinstitute/gatk-bwamem-jni)
[![Build Status](https://travis-ci.org/broadinstitute/gatk-bwamem-jni.svg?branch=master)](https://travis-ci.org/broadinstitute/gatk-bwamem-jni)

# gatk-bwamem-jni
JNI code for bwa mem.

This project builds dynamic libraries with a JNI API.
It allows Java code to call Heng Li's bwa mem aligner.

To build you'll need gmake, git, gcc, and Java 8.

#### To build and install a snapshot locally:

```
./gradlew install
```

This will work for testing but will only include a native library for your system.

#### To upload a snapshot from a Broad Institute OSX machine with both OSX and Linux binaries:
```
commit your changes and push your branch to github
scripts/build_both_dylib_and_so.sh
./gradlew uploadArchives printVersion
```

#### To upload to maven central 

First follow the one-time setup instruction post on the [Picard project Wiki page](https://github.com/broadinstitute/picard/wiki/How-to-release-Picard#one-time-setup-tasks). You can ignore the last part about setting up sudo rights on picard02.

Then do the following:

```
commit your changes and push your branch to github
git tag -a -s <version>
./gradlew clean
scripts/build_both_dylib_and_so.sh
./gradlew uploadArchive -Drelease=true
```
Then you need to finalize the release to maven central:

 * Go to [https://oss.sonatype.org/#stagingRepositories](https://oss.sonatype.org/#stagingRepositories), logging in if necessary. If you don't see anything, click "refresh".

 * Find the release you just uploaded. It will probably be at the bottom with a name like comgithubbroadinstitute-1027, with your user ID listed as owner.

 * Check the box next to your release, then select "close". Press the refresh button repeatedly until it updates and its status says "closed".

 * Select your release again and click "release". Select the box to "automatically drop" in the pop-up confirmation dialog.

It may take ~30-180 minutes for the new release to show up on the maven central.

#### To use this JNI binding on another architecture for which we don't provide a binary:

  Go into ```src/main/c```.
  Modify the Makefile to produce a library name appropriate to your system.
  Type ```make``` (you'll need gmake, git, and gcc).
  Move the library you built somewhere permanent on your machine.
  Use ```-DLIBBWA_PATH=<that permanent location>``` when you run GATK (or other Java program).
