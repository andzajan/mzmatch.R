==============================================
mzMatch/PeakML: metabolomics data analysis 
==============================================

|Git| |Git_commit| |Life_cycle| |License| |mzmatch| |mzMatch-ISO| 

mzMatch is a modular, open source and platform independent data processing pipeline for metabolomics LC/MS data written in the Java language. It was designed to provide small tools for the common processing tasks for LC/MS data. The mzMatch environment was based entirely on the PeakML file format and core library, which provides a common framework for all the tools.

With the introduction of the PeakML file format, the tools can share their data output with the other software. This means that we can integrate with other tools, for which we provide a first implementation integration mzMatch with XCMS. 

This repository contains R package designed to extend functionality of mzmatch and PeakML
modules in the single unified interface.

------------
Install
------------

Prerequisites
------------

- R version 3.6.2 or higher
- Oracle Java 8 JRE (Java Runtime Environment)
- Rstudio (optional)

Installing mzmatch.R package
------------

.. code-block:: r

  install.packages ("remotes")
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  remotes::install_github("https://github.com/andzajan/mzmatch.R.git")

Installing mzmatch java libraries
------------

.. code-block:: r
  library(mzmatch.R)
  mzmatch.init(version.1=FALSE)
  
Installation of PeakML.Viewer
------------

The PeakML.Viewer will be installed when called for the
first time

.. code-block:: r
  PeakML.Viewer()

------------
References
------------

.. |Git| image:: https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600
   :target: https://github.com/computational-metabolomics/qcrms

.. |License| image:: https://img.shields.io/badge/licence-GNU_v3-teal.svg?style=flat&maxAge=3600
   :target: https://www.gnu.org/licenses/gpl-3.0.html

.. |mzMatch| image:: https://img.shields.io/badge/doi-10.1021/ac2000994-yellow.svg
   :target: https://doi.org/10.1021/ac2000994

.. |mzMatch-ISO| image:: https://img.shields.io/badge/doi-10.1093/bioinformatics/bts674-yellow.svg
   :target: https://doi.org/10.1093/bioinformatics/bts674

.. |Git_commit| image:: https://doi.org/10.1021/ac2000994
   :target: https://github.com/andzajan/mzmatch.R/commits/master
   
.. |Life_cycle| image::https://img.shields.io/badge/lifecycle-stable-green.svg
   :target: https://www.tidyverse.org/lifecycle/#stable
   