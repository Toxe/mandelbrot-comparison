all: release

debug:
	swiftc -g main.swift mandelbrot.swift -o mandelbrot

release:
	swiftc -O main.swift mandelbrot.swift -o mandelbrot

clean:
	rm -rf mandelbrot
