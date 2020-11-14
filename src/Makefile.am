ACLOCAL_AMFLAGS= -I m4
# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@

@VALGRIND_CHECK_RULES@

CLEANFILES =

AM_CFLAGS = @TLDEVEL_CPPFLAGS@

AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef
AM_LDFLAGS += -static

LIBS = @LIBS@ @TLDEVEL_LIB@ -lm

bin_PROGRAMS = phenoparser

phenoparser_SOURCES = \
phenoparser.h \
omimlist.c \
parser.c \
sqlite_database.c \
main.c \
user_interface.c \
misc.c
