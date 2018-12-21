#!/bin/csh -f
#
# Shell file for generating pstex files (for the first time)
#
# Usage:
#   see genpstex.sh for usage details
# Created:
#   Chris Bradley
# Updates:
#   Chris Braldey 10/3/96 Using the same file for plots and figs
#
if($1 == 'plots') then
make -f ${OPENCMISS_ROOT}/src/iron/doc/latex/Makeplots $2:r.pstex
else
make -f ${OPENCMISS_ROOT}/src/iron/doc/latex/Makefigs $2:r.pstex
endif
