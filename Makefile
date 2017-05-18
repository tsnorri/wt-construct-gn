include local.mk
include common.mk

DEPENDENCIES = lib/sdsl/build/lib/libsdsl.a
ifeq ($(shell uname -s),Linux)
	DEPENDENCIES    +=  lib/libdispatch/libdispatch-build/src/libdispatch.a
	DEPENDENCIES    +=  lib/libpwq/libpwq-build/libpthread_workqueue.a
endif


.PHONY: all clean-all clean clean-dependencies dependencies


all: dependencies
	$(MAKE) -C src

clean-all: clean clean-dependencies

clean:
	$(MAKE) -C src clean

clean-dependencies:
	cd lib/sdsl/build && ./clean.sh
	$(RM) -r lib/libdispatch/libdispatch-build
	$(RM) -r lib/libpwq/libpwq-build

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

lib/libdispatch/libdispatch-build/src/libdispatch.a: lib/libpwq/libpwq-build/libpthread_workqueue.a
	rm -rf lib/libdispatch/libdispatch-build && \
	cd lib/libdispatch && \
	mkdir libdispatch-build && \
	cd libdispatch-build && \
	../configure --cc="$(CC)" --c++="$(CXX)" --release -- \
		-DPTHREAD_WORKQUEUE_INCLUDE_DIRS=../../libpwq/include \
		-DPTHREAD_WORKQUEUE_LIBRARIES=../../libpwq/libpwq-build/libpthread_workqueue.a \
		-DBLOCKS_RUNTIME_LIBRARIES=""
	$(MAKE) -C lib/libdispatch/libdispatch-build VERBOSE=1

lib/libpwq/libpwq-build/libpthread_workqueue.a:
	rm -rf lib/libpwq/libpwq-build && \
	cd lib/libpwq && \
	mkdir libpwq-build && \
	cd libpwq-build && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	cmake -DSTATIC_WORKQUEUE=ON ..
	$(MAKE) -C lib/libpwq/libpwq-build VERBOSE=1
