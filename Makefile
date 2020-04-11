# -*- mode: makefile; -*-

include $(wildcard config/*.mk)

TARGET = target
BUILD_DIR = build
VALGRIND_DIR = build/codetest
DOXYGEN_DIR = doc
COVERAGE_DIR = build
INSTALL_DIR = $(HOME)/local
ARG =  -Db_coverage=true $(ASAN) -Dprefix=$(INSTALL_DIR)
PKGLIB = 0
BUILD = -DUSE_PKG_LIB=$(PKGLIB) -DNO_COTIRE=1 -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY=ON -Wno-dev 
ifneq ($(NINJA),)
  BUILD := $(BUILD) -GNinja
endif

# R package:
pkg = gof
TESTR = $(pkg)_test

# python package
TESTPY = targeted_test

##################################################

default: checkinit build runr

all: clean run

##################################################

.PHONY: clean
clean: cleanr cleanpy
	@rm -Rf $(BUILD_DIR) $(VALGRIND_DIR) $(DOXYGEN_DIR)/html $(COVERAGE_DIR)
	@rm -Rf src/*.o src/*.so

.PHONY: init init-submodules checkinit
init: clean
	@echo "Build options: $(BUILD)"
	@echo $(shell pwd)
	@mkdir -p build
	@echo "Comfiguration..."
	@cd build; $(CMAKE) .. $(BUILD)

checkinit:
	@if [ ! -f "$(BUILD_DIR)/build.ninja" ]; then $(MAKE) init; fi

init-submodules:
	@if [ -z "`find \"lib/armadillo\" -mindepth 1 -exec echo notempty \; -quit`" ]; then \
	$(GIT) submodule update --init --recursive; fi

.PHONY: run
run: 
	@if [ ! -d "$(BUILD_DIR)" ]; then $(MAKE) --no-print-directory init; fi
	@$(MAKE) --no-print-directory build # > /dev/null
	@printf "\n-----\n"
	@find build/ -maxdepth 1 -iname "*demo" $(FINDEXEC) \
	-exec {} \; 

.PHONY: build
build:
	@if [ -f $(BUILD_DIR)/build.ninja ]; then \
	$(NINJA) -C $(BUILD_DIR) $(NINJA_BUILD_OPT); \
	else \
	cd $(BUILD_DIR); make; \
	fi

.PHONY: install
install:
	$(NINJA) -C $(BUILD_DIR) install

.PHONY: uninstall
uninstall:
	$(NINJA) -C $(BUILD_DIR) uninstall

##################################################
## R package
##################################################

.PHONY: r cleanr buildr runr testr roxygen
buildr: cleanr
	@$(R) --slave -e "source('config/utilities.R'); \
	load_packages(c('Rcpp', 'RcppArmadillo', 'lava', 'optimx', 'futile.logger'))"
	@$(R) --slave -e "Rcpp::compileAttributes('R-package/${pkg}')"
	@$(R) CMD INSTALL R-package/${pkg}

testr:
	@$(R) -e 'testthat::test_package("./R-package/${pkg}/")'

runr:
	@cd misc; $(R) --silent -f $(TESTR).R

roxygen:
	@$(R) -e 'roxygen2::roxygenize("R-package/${pkg}")'

dep_file = R-package/${pkg}/src/dependencies
pkg_dep := $(shell if [ -f "$(dep_file)" ]; then cat ${dep_file}; fi)
pkg_cpp = $(foreach module, ${pkg_dep}, $(patsubst %, src/%.cpp, $(module)))
pkg_hpp = $(foreach module, $(pkg_dep), $(patsubst %, include/target/%.hpp, $(module)))

exportr:
	@rm -Rf $(BUILD_DIR)/R/$(pkg)
	@mkdir -p $(BUILD_DIR)/R/$(pkg)
	cd R-package/${pkg}; $(GIT) archive HEAD | (cd ../../$(BUILD_DIR)/R/$(pkg); tar x)
	@if [ -z "$(pkg_dep)" ]; then \
	cp src/*.cpp $(BUILD_DIR)/R/$(pkg)/src; \
	cp -a include/target/ $(BUILD_DIR)/R/$(pkg)/inst/include/; \
	else \
	cp $(pkg_cpp) $(BUILD_DIR)/R/$(pkg)/src; \
	mkdir -p $(BUILD_DIR)/R/$(pkg)/inst/include/target; \
	cp $(pkg_hpp) $(BUILD_DIR)/R/$(pkg)/inst/include/target; \
	fi
	sed -i '/^OBJECTS\|SOURCES/d' $(BUILD_DIR)/R/$(pkg)/src/Makevars
	cd $(BUILD_DIR)/R; $(R) CMD build $(pkg) --compact-vignettes=gs+qpdf --resave-data=best

checkr: exportr
	cd $(BUILD_DIR)/R; $(R) CMD check `$(GETVER) $(pkg)` --timings --as-cran --no-multiarch --run-donttest

r: buildr runr

cleanr:
	@rm -Rf R-package/${pkg}/src/*.o R-package/${pkg}/src/*.so

.PHONY: syncr
syncr: exportr
	@if [ -d "../$(pkg)" ]; then \
	cp -rfv $(BUILD_DIR)/R/$(pkg)/.	 ../$(pkg)/; fi

##################################################
## Python package
##################################################

.PHONY: py cleanpy buildpy runpy testpy

buildpy:
	@cd python-package; $(PYTHON) setup.py install

testpy:
	@cd python-package; $(MAKE) test

runpy:
	@$(PYTHON) misc/$(TESTPY).py

py: buildpy runpy

PYTHON_EXPORT = $(BUILD_DIR)/python
exportpy: clean
	@rm -Rf $(PYTHON_EXPORT); mkdir -p $(PYTHON_EXPORT)
	@cd python-package; $(GIT) archive HEAD | (cd ../$(PYTHON_EXPORT); tar x)
	@cp -a src $(PYTHON_EXPORT)/lib/target-cpp
	@cp -a include $(PYTHON_EXPORT)/lib/target-inc
	@cp -a lib/armadillo $(PYTHON_EXPORT)/lib
	@cp -a lib/catch2 $(PYTHON_EXPORT)/lib
	@cp -a lib/pybind11 $(PYTHON_EXPORT)/lib
	@echo "Python package exported to: ${PYTHON_EXPORT}"

cleanpy:
	@cd python-package; $(MAKE) --no-print-directory clean

##################################################
## Documentation
##################################################

.PHONY: docs doc
docs:
	@cd $(DOXYGEN_DIR); $(DOXYGEN)

doc:	docs
	@$(OPEN) $(DOXYGEN_DIR)/html/index.html

markdown:
	@if [ -z "command -v grip" ]; then \
	echo "Install dependency: pip install grip"; \
	else grip 1313 -b; fi

##################################################
## Unit tests
##################################################

.PHONY: t test testall
t:	run
	@$(NINJA) -C $(BUILD_DIR) test

test:	build
	build/$(TARGET)_test -s

testall: test r py testr testpy

##################################################
## Code coverage
##################################################
.PHONY: coverage
coverage:
	rm -Rf build/coverage
	mkdir -p build/coverage
	cd build/coverage; cmake -DCMAKE_BUILD_TYPE=Debug -DCOVERAGE_BUILD=1 ../../ && make coverage
	$(OPEN) build/coverage/coverage/index.html

##################################################
## Debugging, profiling, and memory inspection
##################################################

.PHONY: check
check:
	-cclint src/*.cpp src/*.h
	-cppcheck --enable=all src/

.PHONY: valgrind
## Alternatively, enable Address Sanitizer (ASAN =-Db_sanitize=address)
valgrind:
	@$(NINJA) -C build test_memcheck

##################################################
## Container
##################################################

.PHONY: dockerbuild dockerrun docker export
dockerbuild:
	@$(CONTAINER_RUNTIME) build . --network=host -t $(TARGET)

export:
	@rm -Rf ${PWD}/tmp/$(TARGET)
	@mkdir -p ${PWD}/tmp/$(TARGET)
	@git archive HEAD | (cd ${PWD}/tmp/$(TARGET); tar x)
	@git submodule foreach 'curdir=${PWD} cd ${PWD}/$$path; git archive HEAD | tar -x -C ${PWD}/tmp/$(TARGET)/$$path'
	@echo "Exported to '${PWD}/tmp/$(TARGET)'"
	@chmod -R 777 ${PWD}/tmp/$(TARGET)

dockerrun:
	$(CONTAINER_RUNTIME) run --user `id -u` -ti --rm --privileged -v ${PWD}/tmp/$(TARGET):/data $(TARGET) ${CMD}

docker: export dockerrun
