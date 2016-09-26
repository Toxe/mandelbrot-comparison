all: debug

debug: mandelbrot.c
	cc -Wall -O0 -g -o mandelbrot mandelbrot.c

release: mandelbrot.c
	cc -Wall -O3 -o mandelbrot mandelbrot.c

clean:
	rm -r mandelbrot *.dSYM
