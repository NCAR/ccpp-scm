DOCKER ?= docker
CONTEXT ?= ..
NETCDF_FORTRAN_VERSION ?= 4.6.2
COMPILER := gnu
# COMPILER := oneapi
# COMPILER := nvhpc
FINAL_IMAGE ?= scm-$(COMPILER):latest

.PHONY: all build minimal netcdf nceplibs python finalize run clean

all: build

build: finalize dev

minimal:
	$(DOCKER) build \
		-f Dockerfile-$(COMPILER)-minimal \
		-t minimal-toolchain:$(COMPILER) \
		$(CONTEXT)

netcdf: minimal
	$(DOCKER) build \
		-f Dockerfile-add-netcdf \
		--build-arg TOOLCHAIN=$(COMPILER) \
		--build-arg NFVERSION=$(NETCDF_FORTRAN_VERSION) \
		-t add-netcdf:$(COMPILER) \
		$(CONTEXT)

nceplibs: netcdf
	$(DOCKER) build \
		-f Dockerfile-add-nceplibs \
		--build-arg TOOLCHAIN=$(COMPILER) \
		-t add-nceplibs:$(COMPILER) \
		$(CONTEXT)

python: nceplibs
	$(DOCKER) build \
		-f Dockerfile-add-python \
		--build-arg TOOLCHAIN=$(COMPILER) \
		-t add-python:$(COMPILER) \
		$(CONTEXT)

finalize: python
	$(DOCKER) build \
		-f Dockerfile-finalize \
		--build-arg TOOLCHAIN=$(COMPILER) \
		-t $(FINAL_IMAGE) \
		$(CONTEXT)

dev: python
	$(DOCKER) build \
		-f Dockerfile-dev \
		--build-arg TOOLCHAIN=$(COMPILER) \
		-t scm-$(COMPILER):dev \
		$(CONTEXT)

run:
	@echo "From the top directory run one of the following:"
	@echo " - Image for building and running SCM"
	@echo "$$ docker run --rm -it --user "$$(id -u):$$(id -g)" --env HOME=/tmp -v "$${PWD}/../:/src" -w /src scm-$(COMPILER):latest"
	@echo " - Image for SCM development"
	@echo "$$ docker run --rm -it --user "$$(id -u):$$(id -g)" --env HOME=/tmp -v "$${PWD}/../:/src" -w /src scm-$(COMPILER):dev"

clean:
	$(DOCKER) image rm \
		$(FINAL_IMAGE) \
		scm-$(COMPILER):dev \
		add-python:$(COMPILER) \
		add-nceplibs:$(COMPILER) \
		add-netcdf:$(COMPILER) \
		minimal-toolchain:$(COMPILER)
