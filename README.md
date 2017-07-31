bacter
======

[![Build Status](https://travis-ci.org/tgvaughan/bacter.svg?branch=master)](https://travis-ci.org/tgvaughan/bacter)
[![CodeFactor](https://www.codefactor.io/repository/github/tgvaughan/bacter/badge)](https://www.codefactor.io/repository/github/tgvaughan/bacter)

Bacter is a [BEAST 2](http://www.beast2.org)  package which facilitates
inference of a (restricted kind of) ancestral recombination graph (ARG) and
related parameters from a sequence alignment.  It is based on the model
described in [Didelot et al.'s 2010 Genetics paper][1].

This archive contains the source code of the package and is therefore of
primary interest to programmers.  For installation and usage instructions, as
well as links to tutorials and other documentation, please visit the project
home page hosted at http://tgvaughan.github.io/bacter.

Archive Contents
----------------

* `/examples` : Example XML files, simulated data for the tutorial and a
  Jupyter notebook with implementation validation details. (You can view this
  notebook [online][2].) 
* `/lib` : Required libraries.
* `/src` : Java source code.
* `/test` : Java source code (unit tests).
* `/templates` : BEAUti templates.
* `version.xml` : BEAST package version file.
* `build.xml` : Ant build script.
* `.gitignore` : Causes SCM to ignore certain files
* `.travis.yml` : Control file for Travis CI server
* `Dockerfile`: Used by Travis to build a reproducible test environment.
* `COPYING` : Software license.
* `README.md` : This file.

Building package from source
----------------------------

To build this package from source, ensure you have the following installed:

* Java JDK v1.8 
* Apache Ant v1.9 or later
* An internet connection

The internet connection is required since the build script downloads the most
recent version of the BEAST 2 source to build the package against.
Assuming both Java and Ant are on your execution path and your CWD is the root of
this archive, simply type "ant" from the command line to build the package.
This may take up to a minute due to the script fetching the BEAST source, and
the resulting binary will be left in the `/dist` directory.
To run the unit tests, use "ant test".

License
-------

This software is free (as in freedom). You are welcome to use it, modify it,
and distribute your modified versions provided you extend the same courtesy to
users of your modified version.  Specifically, it is made available under the
terms of the GNU General Public License version 3, which is contained in his
directory in the file named COPYING.

Acknowledgements
----------------

Work on this project is made possible by the support of the following institutions:

* [The Allan Wilson Centre for Molecular Ecology and Evolution](http://www.allanwilsoncentre.ac.nz)

* [The Royal Society of New Zealand's Marsden Fund](http://www.royalsociety.org.nz/programmes/funds/marsden/)

* [The University of Auckland](http://auckland.ac.nz)

* [Massey University](http://www.massey.ac.nz)

[1]: http://www.genetics.org/content/186/4/1435
[2]: http://nbviewer.jupyter.org/github/tgvaughan/bacter/blob/master/examples/Validation.ipynb
