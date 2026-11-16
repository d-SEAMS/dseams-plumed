# Same contract as PLUMED's regtest makefiles: `make test` in an rt-* directory.
ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../..)
.PHONY: test
test:
	@bash $(ROOT)/regtest/run
