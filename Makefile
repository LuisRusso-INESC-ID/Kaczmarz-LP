# lpgraph build.  DEPS points at the directory where HiGHS, GLPK, Clp (and optionally
# OR-Tools under $(DEPS)/ortools) are installed; LEMON headers are needed for the baselines.
#
# Two LP binaries are built so that the large OR-Tools shared libraries (abseil, protobuf,
# ~10 MB resident) do not inflate the memory baseline of the other backends:
#   bin/lpsolve      HiGHS, Clp, GLPK and the own prototypes
#   bin/lpsolve-ort  OR-Tools GLOP and PDLP (and the own prototypes)
DEPS   ?= /home/claude/deps/install
LEMON  ?= /home/claude/deps/lemon/include
CXX    ?= g++
CXXFLAGS ?= -O2 -std=c++17 -DNDEBUG
ORT    = $(DEPS)/ortools

all: bin/gen bin/classic bin/lpsolve $(if $(wildcard $(ORT)/include/ortools/pdlp),bin/lpsolve-ort)

HDRS = src/lp.h src/pdhg.h src/skm.h src/prowact.h src/hybrid.h src/bisect.h src/instance.h

bin/gen: src/gen.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
bin/classic: src/classic.cpp src/instance.h
	$(CXX) $(CXXFLAGS) -I$(LEMON) -o $@ $<
bin/lpsolve: src/lpsolve.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -march=native -DHAVE_HIGHS -DHAVE_GLPK -DHAVE_CLP \
	  -I$(DEPS)/include/highs -I$(DEPS)/include -I$(DEPS)/include/coin-or -Isrc -o $@ $< \
	  -L$(DEPS)/lib -Wl,-rpath,$(DEPS)/lib -lhighs -lglpk -lClp -lCoinUtils
# no -march=native here: Eigen objects cross the OR-Tools ABI boundary
bin/lpsolve-ort: src/lpsolve.cpp src/ortools_backend.h $(HDRS)
	$(CXX) $(CXXFLAGS) -DHAVE_ORTOOLS -DOR_PROTO_DLL= -DUSE_PDLP -DUSE_GLOP \
	  -I$(ORT)/include -I$(ORT)/include/eigen3 -Isrc -o $@ $< \
	  -L$(ORT)/lib -Wl,-rpath,$(ORT)/lib -lortools -lprotobuf \
	  $(patsubst $(ORT)/lib/lib%.so,-l%,$(wildcard $(ORT)/lib/libabsl_*.so))

clean:
	rm -f bin/*
.PHONY: all clean
