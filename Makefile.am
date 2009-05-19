ACLOCAL_AMFLAGS=-I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libspace.pc

etags:
	find . -name "*.cc" -or -name "*.h" -or -name "*.c" -or -name "*.cpp" -or -name "*.C" | xargs etags -l c++


#clean:
#	find . -\( -name "*~" -or -name "\#*\#" -or -name "*.o" -or -name ".deps" -or -name ".libs" -or -name "*.Plo" -or -name "*.lo" -\) -exec rm -rf {} \;
#
clean-output:
	rm -f *.vtk *.gpl *.eps


SUBDIRS=src test
