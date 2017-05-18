include local.mk
include common.mk

DEPENDENCIES = lib/sdsl/build/lib/libsdsl.a

all: dependencies
	$(MAKE) -C src

clean-all: clean clean-dependencies

clean:
	$(MAKE) -C src clean

clean-dependencies:
	cd lib/sdsl/build && ./clean.sh

dependencies: $(DEPENDENCIES)

lib/sdsl/build/lib/libsdsl.a:
	cd lib/sdsl/build && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	CFLAGS="$(CFLAGS) $(CPPFLAGS) $(OPT_FLAGS)" \
	CXXFLAGS="$(CXXFLAGS) $(CPPFLAGS) $(OPT_FLAGS)" \
	LDFLAGS="$(LDFLAGS) $(OPT_FLAGS)" \
	cmake ..
	$(MAKE) -C lib/sdsl/build VERBOSE=1
