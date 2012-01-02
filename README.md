Overview
========

**Genherite** predicts genetic networks for multicellular cell fate decisions.

A distribution version is available as a
[tarball](https://github.com/hrouault/Genherite/downloads).

It is programmed in C/C++. Source files can be found in the folder `src`.


Installation
============

You can refer to the `INSTALL` file for detailled installation instructions. This 
file is readily available in the tarball. If you use the development version, you will
need the autotools in order to generate it and compile the code. 

Quick install
-------------

```sh
mkdir build
cd build
../configure
make
make install
```

You can control the destination folder by adjusting the `prefix` variable at
the configure step:

```sh
../configure --prefix=/folder/where/to/install
```

see the `INSTALL`file for more informations.

Development version
-------------------

Run:

```sh
./autogen.sh
```

prior to the steps described in the previous paragraph.

Running Genherite
=================

Run `src/genherite --help` to see the different modes of execution.
