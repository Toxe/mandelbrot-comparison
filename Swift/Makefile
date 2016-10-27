all: release

debug:
	swiftc -g main.swift mandelbrot.swift -o mandelbrot

release:
	swiftc -O -whole-module-optimization main.swift mandelbrot.swift -o mandelbrot

clean:
	rm -rf mandelbrot
