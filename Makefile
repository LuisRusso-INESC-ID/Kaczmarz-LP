# lpgraph build.  DEPS points at the directory where HiGHS, GLPK, Clp (and optionally
# OR-Tools under $(DEPS)/ortools) are installed; LEMON headers are needed for the baselines.
#
# Two LP binaries are built so that the large OR-Tools shared libraries (abseil, protobuf,
# ~10 MB resident) do not inflate the memory baseline of the other backends:
#   bin/lpsolve      HiGHS, Clp, GLPK and the own prototypes
#   bin/lpsolve-ort  OR-Tools GLOP and PDLP (and the own prototypes)
DEPS   ?= $(HOME)/lpdeps/install
LEMON  ?= $(HOME)/lpdeps/src/lemon/include
CXX    ?= g++
CXXFLAGS ?= -O2 -std=c++17 -DNDEBUG
ORT    = $(DEPS)/ortools
# reference max-flow codes, fetched by scripts/fetch_flow_codes.sh (not redistributed; see there)
FLOW   ?= $(HOME)/lpdeps/src/flow
FLOWBIN = $(if $(wildcard $(FLOW)/hipr42/main.c),bin/flow-hipr) \
          $(if $(wildcard $(FLOW)/pseudo_max/src/3.23/pseudo.c),bin/flow-pseudoflow) \
          $(if $(wildcard $(FLOW)/maxflow-v3.04.src/maxflow.cpp),bin/flow-bk)

all: bin/gen bin/classic bin/lpsolve $(if $(wildcard $(ORT)/include/ortools/pdlp),bin/lpsolve-ort) $(FLOWBIN)

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

# authors' optimization flags; -std=gnu17 (GCC 15+ defaults to C23) and -include stdlib.h
# and string.h (the 2008-era C calls malloc/free/strcmp undeclared, an error in GCC 14+)
bin/flow-hipr: src/flowref/hipr_driver.c
	$(CC) -O3 -std=gnu17 -include stdlib.h -include string.h -DNDEBUG -DFULL_RELABEL -DLASY_UPDATE -w -I$(FLOW)/hipr42 -o $@ $< \
	  $(FLOW)/hipr42/parser.c $(FLOW)/hipr42/timer.c -lm
bin/flow-pseudoflow: src/flowref/pseudo_driver.c
	$(CC) -O3 -std=gnu17 -include stdlib.h -DLOWEST_LABEL -DFIFO_BUCKET -w -I$(FLOW)/pseudo_max/src/3.23 -o $@ $<
bin/flow-bk: src/flowref/bk_driver.cpp src/flowref/bk_maxflow.cpp
	$(CXX) -O3 -std=c++17 -DNDEBUG -w -I$(FLOW)/maxflow-v3.04.src -o $@ $^

clean:
	rm -f bin/*
.PHONY: all clean
