#-----------
# make-include.mak--
#	Defines variables used by all LASTZ Makefiles
#-----------

INSTALL =  install
ARCH    ?= $(shell uname -m)

ifdef LASTZ_INSTALL
installDir = ${LASTZ_INSTALL}
else
installDir = ${HOME}/lastz-distrib/bin
endif

