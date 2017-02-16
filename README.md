# gatk-bwamem-jni
JNI code for bwa mem.

This project builds dynamic libraries with a JNI API.
It allows Java code to call Heng Li's bwa mem aligner.

The makefile in src/main/c will build an appropriate library for Mac OSX or x86_64 Linux.
Pre-compiled dynamic libraries for these OS's exist in src/main/resources.

To deploy into maven central:
Find a Mac.
  Clone this repository.
  Go into ```src/main/c```.
  Type ```make``` (you'll need gmake, git, and gcc).
  Copy ```libbwa.Darwin.dylib``` to ```src/main/resources```.
  Type ```make clean```.
  Go back up to the repo root, and type ```gradle test``` to run the Java unit tests.
Find a Linux machine.
  Clone this repository.
  Go into ```src/main/c```.
  Type ```make``` (you'll need gmake, git, and gcc).
  Copy ```libbwa.Linux.so``` to ```src/main/resources```.
  Type ```make clean```.
  Go back up to the repo root, and type ```gradle test``` to run the Java unit tests.
  Copy ```src/main/resources/libbwa.Linux.so``` to ```src/main/resources``` ON THE MAC.
    (Yes, you have to copy it to the other machine.)
Go back to the Mac.
  Type ```gradle deploy``` (you'll need gradle).


To use this JNI binding on some architecture for which we don't provide a binary:
  Clone the repo.
  Go into ```src/main/c```.
  Modify the Makefile to produce a library name appropriate to your system.
  Type ```make``` (you'll need gmake, git, and gcc).
  Move the library you built somewhere permanent on your machine.
  Use ```-DLIBBWA_PATH=<that permanent location>``` when you run GATK (or other Java program).

