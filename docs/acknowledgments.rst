.. _acknowledgments:

###############
Acknowledgments
###############

*******
Funding
*******

This software benefited from funding provided to Jamie Oaks from the National
Science Foundation (DBI 1308885 and DEB 1656004).

****
Code
****

|Eco| comes bundled with the `Nexus Class Library
<http://ncl.sourceforge.net/>`_
Version 2.1.18 (Copyright Paul O. Lewis and Mark T. Holder; licensed under a
BSD License).
NCL is used for parsing the nexus input filesa.

|Eco| also comes bundled with
`yaml-cpp <https://github.com/jbeder/yaml-cpp>`_
commit 7d2873c (master branch) (Copyright Jesse Beder; licensed under a
BSD-like License).
Yaml-cpp is used for parsing |eco| configuration files.

|Eco| also comes bundled with
`cpp-optparse <https://github.com/weisslj/cpp-optparse>`_
(Copyright Johannes Weiß; MIT License).
Cpp-optparse is used to parse command-line options.

The likelihood machinery in |eco| was adapted from the java code in
`SNAPP <http://www.beast2.org/snapp/>`_
which is a package for
`BEAST 2 <http://www.beast2.org/>`_
(GNU Lesser General Public License Version 2.1).
Citation: Bryant, D., Bouckaert, R., Felsenstein, J., Rosenberg, N. A., and
Roychoudhury, A. 2012.  Inferring species trees directly from biallelic genetic
markers: Bypassing gene trees in a full coalescent analysis.  Molecular Biology
and Evolution, 29(8): 1917–1932.  https://doi.org/10.1093/molbev/mss086

Some of the convenience functions for calculating properties of Dirichlet
processes were inspired by functions within ``util.h`` from the software
package `DPPDiv <http://phylo.bio.ku.edu/content/tracy-heath-dppdiv>`_ Version
1.0b (Copyright Tracy Heath, Mark Holder, and John Huelsenback; licensed under
GPL v3; http://phylo.bio.ku.edu/content/tracy-heath-dppdiv).
