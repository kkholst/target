# -*- mode: makefile; -*-

include $(wildcard config/*.mk)
TARGET = target
BUILD_DIR = build
VALGRIND_DIR = build/codetest
INSTALL_DIR = $(HOME)/local
PKGLIB = ON
IMG=# Dockerfile postfix
CXX=
EXTRA=
BUILD=-DUSE_PKG_LIB=$(PKGLIB) -DCOTIRE=OFF -DCMAKE_INSTALL_PREFIX=$(INSTALL_DIR) -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY=ON -Wno-dev -DCMAKE_VERBOSE_MAKEFILE=ON \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=YES -H.
BUILD+=$(EXTRA)
ifneq ($(CXX),)
   BUILD+= -DCMAKE_CXX_COMPILER=$(CXX)
endif
ifneq ($(NINJA),)
   BUILD+= -GNinja
endif

# R package:
pkg = targeted
TESTR = $(pkg)_test

# python package
TESTPY = $(pkg)_test

##################################################

default: build run

all: clean run

##################################################

.PHONY: cleansrc
cleansrc:
	@rm -Rf $(BUILD_DIR) $(VALGRIND_DIR)
	@rm -Rf src/*.o src/*.so

.PHONY: clean
clean: cleansrc

.PHONY: cleanall
cleanall: cleansrc py_clean cleandoc
	@$(MAKE) pkg=gof r_clean
	@$(MAKE) pkg=targeted r_clean

.PHONY: init
init: cleansrc
	@echo "Build options: $(BUILD)"
#	@$(CMAKE) $(BUILD) -B$(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR); $(CMAKE) $(BUILD) ..
	@ln -fs $(BUILD_DIR)/compile_commands.json

.PHONY: checkinit
checkinit: init-submodules
	@if [ ! -f "$(BUILD_DIR)/CMakeCache.txt" ]; then $(MAKE) init; fi

.PHONY: init-submodules
init-submodules:
	@if [ ! -e "lib/armadillo/.git" ] || [ ! -e "lib/doctest/.git" ] \
	|| [ ! -e "lib/spdlog/.git" ] || [ ! -e "lib/pybind11/.git" ]; then \
	$(GIT) submodule update --init --recursive; fi

.PHONY: run
run:
	@if [ ! -d "$(BUILD_DIR)" ]; then $(MAKE) --no-print-directory init; fi
	@$(MAKE) --no-print-directory build # > /dev/null
		@find $(BUILD_DIR)/ -maxdepth 1 \
	\( -iname "*demo" -o -iname "*run" \) $(FINDEXEC) \
	-exec sh -c "printf '___ {} ____________________________________________________________________\n'; {}" \;

.PHONY: build
build: checkinit
	@if [ -f $(BUILD_DIR)/build.ninja ]; then \
	$(NINJA) -C $(BUILD_DIR) $(NINJA_BUILD_OPT); \
	else \
	cd $(BUILD_DIR); make; \
	fi

.PHONY: install
install: checkinit
	$(NINJA) -C $(BUILD_DIR) install

.PHONY: uninstall
uninstall:
	$(NINJA) -C $(BUILD_DIR) uninstall

##################################################
## R package
##################################################

.PHONY: r_build
r_build: r_clean
	@$(R) --slave -e "source('config/utilities.R'); \
	load_packages(c('Rcpp', 'RcppArmadillo', 'lava', 'optimx', 'futile.logger'))"
	@$(R) --slave -e "Rcpp::compileAttributes('R-package/${pkg}')"
	@$(R) CMD INSTALL R-package/${pkg}

.PHONY: r_test
r_test:
	@$(R) -e 'testthat::test_package("./R-package/${pkg}/")'

.PHONY: r_run
r_run:
	@cd misc; $(R) --silent -f $(TESTR).R

.PHONY: r_doc roxygen
r_doc:
	@$(R) -e 'roxygen2::roxygenize("R-package/${pkg}")'
roxygen: r_doc

dep_file = R-package/${pkg}/src/dependencies
pkg_dep := $(shell if [ -f "$(dep_file)" ]; then cat ${dep_file}; fi)
pkg_cpp = $(foreach module, ${pkg_dep}, $(patsubst %, src/%.cpp, $(module)))
pkg_hpp = $(foreach module, $(pkg_dep), $(patsubst %, include/target/%.hpp, $(module)))

.PHONY: r_export
r_export:
	@rm -Rf $(BUILD_DIR)/R/$(pkg)
	@mkdir -p $(BUILD_DIR)/R/$(pkg)
	cd R-package/${pkg}; $(GIT) archive HEAD | (cd ../../$(BUILD_DIR)/R/$(pkg); tar x )
	@if [ -z "$(pkg_dep)" ]; then \
	cp src/*.cpp $(BUILD_DIR)/R/$(pkg)/src; \
	cp -a include/target $(BUILD_DIR)/R/$(pkg)/inst/include/; \
	else \
	cp $(pkg_cpp) $(BUILD_DIR)/R/$(pkg)/src; \
	mkdir -p $(BUILD_DIR)/R/$(pkg)/inst/include/target; \
	cp $(pkg_hpp) $(BUILD_DIR)/R/$(pkg)/inst/include/target; \
	fi
	sed -i $(SED_NOBACKUP) '/^OBJECTS/d' $(BUILD_DIR)/R/$(pkg)/src/Makevars
	sed -i $(SED_NOBACKUP) '/^SOURCES/d' $(BUILD_DIR)/R/$(pkg)/src/Makevars

.PHONY: r_crancheck
r_crancheck: r_export
	@$(R) --slave -e "source('config/utilities.R'); \
	load_packages(c('devtools'))"
	@$(R) -e "pkgbuild::build('$(BUILD_DIR)/R/$(pkg)', args='--compact-vignettes=qpdf --resave-data=best')"
	cd $(BUILD_DIR)/R; $(R) CMD check `$(GETVER) $(pkg)` --timings --as-cran --no-multiarch --run-donttest

.PHONY: r_check
r_check:
	@$(R) --slave -e "Rcpp::compileAttributes('R-package/${pkg}')"
	@cd R-package; $(R) CMD check $(pkg) --no-multiarch

.PHONY: r
r: r_build r_run

.PHONY: r_clean
r_clean:
	@rm -Rf R-package/${pkg}/src/*.o R-package/${pkg}/src/*.so R-package/${pkg}.Rcheck

.PHONY: r_sync
r_sync: r_export
	@if [ -d "../$(pkg)" ]; then \
	cp -rfv $(BUILD_DIR)/R/$(pkg)/.	 ../$(pkg)/; fi

##################################################
## Python package
##################################################

.PHONY: py_build
py_build:
	@cd python-package/$(pkg); $(PYTHON) setup.py install

.PHONY: py_test
py_test:
	@cd python-package/$(pkg); $(MAKE) test

.PHONY: py_run
py_run:
	@$(PYTHON) misc/$(TESTPY).py

.PHONY: py
py: py_build py_run

PYTHON_EXPORT = $(BUILD_DIR)/python/$(pkg)
PYTHON_EXPORT_LIBS = $(foreach file, target-cpp target-inc doctest armadillo pybind11, \
	$(patsubst %, $(PYTHON_EXPORT)/lib/%, $(file)))
.PHONY: py_export
py_export: cleansrc py_clean
	@rm -Rf $(PYTHON_EXPORT); mkdir -p $(PYTHON_EXPORT)
	@cd python-package/$(pkg); $(GIT) archive HEAD | (cd ../../$(PYTHON_EXPORT); tar x)
	@rm -Rf $(PYTHON_EXPORT_LIBS)
	@cp -a src $(PYTHON_EXPORT)/lib/target-cpp
	@cp -a include $(PYTHON_EXPORT)/lib/target-inc
	@cp -a lib/armadillo $(PYTHON_EXPORT)/lib
	@cp -a lib/doctest $(PYTHON_EXPORT)/lib
	@cp -a lib/pybind11 $(PYTHON_EXPORT)/lib
	@echo "Python package exported to: ${PYTHON_EXPORT}"

.PHONY: py_clean
py_clean:
	@cd python-package/$(pkg); $(MAKE) --no-print-directory clean

##################################################
## Documentation
##################################################

.PHONY: doc
doc:
	@cd doc; $(MAKE) html
	@cd doc/latex; $(MAKE)

.PHONY: html
html:
	@cd doc; $(MAKE) html

.PHONY: docs
docs:
	sphinx-autobuild --open-browser doc/source doc/build

.PHONY: cleandoc
cleandoc:
	@cd doc; $(MAKE) clean
	@rm -Rf doc/latex doc/html doc/xml

.PHONY: markdown
markdown:
	@if [ -z "command -v grip" ]; then \
	echo "Install dependency: pip install grip"; \
	else grip 1313 -b; fi

##################################################
## Unit tests
##################################################

.PHONY: t
t:	checkinit run
	@$(NINJA) -C $(BUILD_DIR) test

.PHONY: test
test:	checkinit build
	$(BUILD_DIR)/$(TARGET)_test -s

.PHONY: testall
testall: test r py r_test py_test

##################################################
## Code coverage
##################################################

.PHONY: cov
cov: coverage
	cd $(BUILD_DIR)/coverage; $(MAKE) coverage
	$(OPEN) $(BUILD_DIR)/coverage/coverage/index.html

.PHONY: coverage
coverage:
	rm -Rf $(BUILD_DIR)/coverage; mkdir -p $(BUILD_DIR)/coverage
	cd $(BUILD_DIR)/coverage; $(CMAKE) -Wno-dev -DCMAKE_BUILD_TYPE=Debug -DCODE_COVERAGE=ON ../../ && $(MAKE)

##################################################
## Debugging, profiling, and memory inspection
##################################################

.PHONY: check
check:
	-@cclint src/*.cpp include/target/*.h* # misc/*demo.cpp misc/*run.cpp misc/*test.cpp
	-@cppcheck --enable=warning,style,performance,portability,information,missingInclude --language=c++ --std=c11 -Isrc -Iinclude src/

.PHONY: valgrind
valgrind:
	@$(NINJA) -C build test_memcheck

.PHONY: sanitizer
sanitizer:
	$(MAKE) run test EXTRA="-DCMAKE_CXX_FLAGS='-fsanitize=address -fsanitize=leak -fsanitize=undefined -fsanitize=bool,signed-integer-overflow,null,alignment,float-divide-by-zero,bool' -DCMAKE_CXX_COMPILER=clang++"
	./build/sanitizer_check


##################################################
## Container
##################################################

DOCKER=Dockerfile
DOCKERTAG=$(TARGET)
ifdef ($IMG)
	DOCKER := '$(DOCKER).$(IMG)'
	DOCKERTAG := '$(DOCKERTAG).$(IMG)'
endif

.PHONY: dockerbuild dockerrun docker export dockersan
dockerbuild: export
	@$(CONTAINER_RUNTIME) build . -f Dockerfile --network=host -t $(DOCKERTAG)

export:
	rm -Rf $(BUILD_DIR)/export
	$(GIT) clone . $(BUILD_DIR)/export;
	rm -Rf $(BUILD_DIR)/export/lib
	$(GIT) submodule foreach '$(GIT) clone . ../../$(BUILD_DIR)/export/$$path/'

docker:
	$(CONTAINER_RUNTIME) run --user `id -u` -ti --rm --privileged -v ${PWD}:/data $(TARGET) ${CMD}

## Use 'RD' instead of 'R'
dockersan: r_export
	$(CONTAINER_RUNTIME) run -ti --rm --privileged -v ./build/R/:/data \
	rocker/r-devel-san ${CMD}
