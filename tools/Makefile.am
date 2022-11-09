## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

DEPS = $(top_builddir)/src/libafmm.la
LDADD =	$(top_builddir)/src/.libs/libafmm.a

EXTRA_DIST = \
	afmm-tests.c.in \
	afmm-direct-calc.c.in \
	afmm-tree-calc.c.in

bin_PROGRAMS = \
	afmm-tests \
	afmm-tests-float \
	afmm-direct-calc \
	afmm-direct-calc-float \
	afmm-tree-calc \
	afmm-tree-calc-float

BUILT_SOURCES = \
	afmm-tests.c \
	afmm-tests-float.c \
	afmm-direct-calc.c \
	afmm-direct-calc-float.c \
	afmm-tree-calc.c \
	afmm-tree-calc-float.c

CLEANFILES = \
	afmm-tests.c afmm-tests-float.c \
	afmm-direct-calc.c afmm-direct-calc-float.c \
	afmm-tree-calc.c afmm-tree-calc-float.c

afmm-tests.c: afmm-tests.c.in
	./make-source afmm-tests

afmm-tests-float.c: afmm-tests.c.in
	./make-source-float afmm-tests

afmm-tree-calc.c: afmm-tree-calc.c.in
	./make-source afmm-tree-calc

afmm-direct-calc.c: afmm-direct-calc.c.in
	./make-source afmm-direct-calc

afmm-tree-calc-float.c: afmm-tree-calc.c.in
	./make-source-float afmm-tree-calc

afmm-direct-calc-float.c: afmm-direct-calc.c.in
	./make-source-float afmm-direct-calc

clean:
	rm -f $(CLEANFILES)
