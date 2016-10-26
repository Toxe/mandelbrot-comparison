all: release

debug:
	$(MAKE) -C C debug
	$(MAKE) -C Swift debug

release:
	$(MAKE) -C C release
	$(MAKE) -C Swift release

clean:
	$(MAKE) -C C clean
	$(MAKE) -C Swift clean
