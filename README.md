# gatk-bwamem-jni
JNI code for bwa mem.

This project builds dynamic libraries with a JNI API.
It allows Java code to call Heng Li's bwa mem aligner.

To build you'll need gmake, git, gcc, and Java 8.

To build and install a snapshot locally:

```
./gradlew install
```

This will work for testing but will only include a native library for your system.

To upload a snapshot from a Broad Institute OSX machine with both OSX and Linux binaries:
```
commit your changes and push your branch to github
scripts/build_both_dylib_and_so.sh
./gradlew uploadArchives printVersion
```

To upload to maven central
```
commit your changes and push your branch to github
git tag -a -s <version>
scripts/build_both_dylib_and_so.sh
./gradlew uploadArchive -Drelease=true
```

To use this JNI binding on another architecture for which we don't provide a binary:
  Go into ```src/main/c```.
  Modify the Makefile to produce a library name appropriate to your system.
  Type ```make``` (you'll need gmake, git, and gcc).
  Move the library you built somewhere permanent on your machine.
  Use ```-DLIBBWA_PATH=<that permanent location>``` when you run GATK (or other Java program).
