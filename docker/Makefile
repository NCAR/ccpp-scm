DOCKER ?= docker
CONTEXT ?= ..
NETCDF_FORTRAN_VERSION ?= 4.6.2
COMPILER := gnu
# COMPILER := oneapi
# COMPILER := nvhpc
FINAL_IMAGE ?= scm-$(COMPILER):latest
RELEASE_VERSION ?= v8.0.0-rc.1
SOURCE_REVISION ?= $(shell git -C $(CONTEXT) rev-parse HEAD)
RELEASE_IMAGE = dtcenter/ccpp-scm:$(RELEASE_VERSION)

.PHONY: all build minimal netcdf nceplibs python finalize release-check release run clean

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

release-check:
	@test "$$(git -C $(CONTEXT) rev-parse HEAD)" = "$(SOURCE_REVISION)" || { \
		printf '%s\n' 'SOURCE_REVISION must match the checked-out commit.' >&2; \
		exit 1; \
	}
	@test -z "$$(git -C $(CONTEXT) status --porcelain)" || { \
		printf '%s\n' 'Release images must be built from a clean checkout.' >&2; \
		git -C $(CONTEXT) status --short >&2; \
		exit 1; \
	}

release: release-check python
	$(DOCKER) build \
		-f Dockerfile-release \
		--build-arg TOOLCHAIN=$(COMPILER) \
		--build-arg RELEASE_VERSION=$(RELEASE_VERSION) \
		--build-arg SOURCE_REVISION=$(SOURCE_REVISION) \
		-t $(RELEASE_IMAGE) \
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
	@echo " - Image for SCM Release"
	@echo "$$ docker run --rm -it $(RELEASE_IMAGE)"

clean:
	@for img in \
		$(FINAL_IMAGE) \
		$(RELEASE_IMAGE) \
		scm-$(COMPILER):dev \
		add-python:$(COMPILER) \
		add-nceplibs:$(COMPILER) \
		add-netcdf:$(COMPILER) \
		minimal-toolchain:$(COMPILER); \
	do \
		if $(DOCKER) image inspect $$img > /dev/null 2>&1; then \
			$(DOCKER) image rm $$img; \
		fi; \
	done
