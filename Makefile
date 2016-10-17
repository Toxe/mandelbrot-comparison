all: release

debug:
	$(MAKE) -C C debug

release:
	$(MAKE) -C C release

clean:
	$(MAKE) -C C clean
