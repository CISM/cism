if EXT_BLAS
  BLAS = $(BLAS_LIBS)
else
  BLAS = ../../blas/libglimmer-blas.a
endif

lib_LIBRARIES = libslap.a

libslap_a_SOURCES = dir.f dcg.f dcgn.f dbcg.f dcgs.f domn.f dgmres.f dmset.f dmvops.f dlaputil.f xersla.f mach.f
check_PROGRAMS = dlapqc
TESTS = runTests.sh

dlapqc_SOURCES=dlapqc.f
dlapqc_LDADD=libslap.a $(BLAS)
