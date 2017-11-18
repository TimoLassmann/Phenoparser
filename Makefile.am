EXTRA_DIST = makeTAGS.sh autogen.sh
# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 

# top_builddir = @top_builddir@

SUBDIRS = tldevel src


TAGS:
	etags -o TAGS -a ${srcdir}/src/*.c
	etags -o TAGS -a ${srcdir}/src/*.h
	etags -o TAGS -a ${srcdir}/tldevel/*.h
	etags -o TAGS -a ${srcdir}/tldevel/*.c
