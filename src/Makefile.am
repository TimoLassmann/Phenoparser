# Build the library in the hand subdirectory first.


LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 

MYLIBDIRS = ../tldevel/libtldevel.a 


MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
-I${top_srcdir}/${LIB_TLDEVELDIR}

AM_CPPFLAGS=${MYINCDIRS} 

bin_PROGRAMS = phenoparser
phenoparser_SOURCES = \
phenoparser.h \
main.c \
user_interface.c \
misc.c


phenoparser_LDADD = ${MYLIBDIRS}

