all: debug

debug: mandelbrot.c
	cc -Wall -Wfloat-equal -O0 -g -o mandelbrot mandelbrot.c

release: mandelbrot.c
	cc -Wall -Wfloat-equal -O3 -o mandelbrot mandelbrot.c

clean:
	rm -r mandelbrot *.dSYM
