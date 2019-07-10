all: release

debug:
	swiftc -g main.swift mandelbrot.swift -o mandelbrot_swift

release:
	swiftc -O -whole-module-optimization main.swift mandelbrot.swift -o mandelbrot_swift

clean:
	rm mandelbrot_swift
