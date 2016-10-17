# Mandelbrot Comparison
A comparison of Mandelbrot Set programs in different languages with histogram based smooth coloring and a built-in benchmark mode.

## Status

| Language | Status |
| -------- | ------ |
| C        | Done.  |
| Python   | WIP    |
| Swift    | TODO   |
| Java     | Maybe? |

## Some background

I made these programs while learning Swift: I needed a not too simple yet on the other hand not too complex example program to translate into Swift. And I wanted something for performance measurements, to compare Swift with C. That's how this multi-language comparison was born.

## Some thoughts what these programs are about

I want to compare the shape, size and performance of a simple yet computationally intensive program between different languages. The programs should be real life examples without too much fancy stuff.

1. The whole code in one file.
2. No dependencies. Just the language and standard library.
3. Simple code, nothing fancy. That means no assembly optimizations or calling C functions from within Python that do all the work.
4. No multi-threading or parallel code. (At least not for the moment.)
5. Error handling is pretty poor. Basically if something went wrong the program quits with some error code. Although for the purpose of these programs this is fine.
6. Keep memory allocations or file access outside the measurement loop.

## Building

In general simply run `make` or `make release` (or `make debug` for a debug build) from the top directory or the C or Swift subdirectories.

## Running

Run the compiled programs from the command line with appropriate parameters.

Command line arguments (see below for examples):

```
$ mandelbrot <image_width> <image_height> <iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
```

- *image_width*, *image_height*: Width and height in pixels of the output image, for example 800 and 600.
- *iterations*: Maximum number of iterations while calculating the Mandelbrot Set, for example 100.
- *repetitions*: Number of times the whole image should be generated. Pick `1` to simply render one image. Pick a bigger value (for example 10) to enter benchmark mode and repeat and measure the image calculations.
- *center x*, *center y*, *section height*: The center of the image in the Mandelbrot Set coordinate system. The height determines the height of the "window". Pick -0.8, 0.0 and 2.2 for the basic, unzoomed Mandelbrot Set image.
- *gradient filename*: Name of the gradient file that should be read in for coloring the image, for example `grey.gradient`.
- *output filename*: The name of the file where the raw RGB output image data should be saved, for example `mandelbrot.raw`.

### C
```
$ ./C/mandelbrot 800 600 1000 10 -0.8 0.0 2.2 grey.gradient mandelbrot.raw
```
### Swift
```
$ ./Swift/mandelbrot 800 600 1000 10 -0.8 0.0 2.2 grey.gradient mandelbrot.raw
```
### Python
```
$ python ./Python/mandelbrot.py 800 600 1000 10 -0.8 0.0 2.2 grey.gradient mandelbrot.raw
```

## Output

Once the image data has been generated it will be saved under the name of the output file (for example `mandelbrot.raw`). The format is raw 8 bit RGB data.

Also it will print out how long it took to calculate and colorize the image.

### Converting the raw image with ImageMagick

To convert the raw image into a more portable format simply use ImageMagick (or some other program that can read raw RGB data).

```
$ convert -size 800x600 -depth 8 rgb:mandelbrot.raw mandelbrot.png
```

Please note that the `-size` parameter needs to match the *width* and *height* arguments above.

It's easiest to simply combine both calls into one command line:

```
$ ./mandelbrot 800 600 1000 1 -0.8 0.0 2.2 grey.gradient mandelbrot.raw && convert -size 800x600 -depth 8 rgb:mandelbrot.raw mandelbrot.png
```

### Benchmark mode

If the *repetitions* command line argument is bigger than one the program will enter benchmark mode. This simply means that it will repeat calculating and colorizing the image as much times as specified in *repetitions*. At the end it will print out the mean and median times and a sorted list of all measurements.

```
mean: 0.897925 s, median: 0.897584 s (repetitions=10) [0.893949, 0.895080, 0.896350, 0.896801, 0.896827, 0.898341, 0.898766, 0.899246, 0.899974, 0.903912]
```

## Gradient files

Files used to colorize the image. Every point in the image gets mapped to a value from 0.0 to 1.0 depending on the number of iterations it took to bail out of the Mandelbrot calculation loop. This value defines a position in a color gradient.

Example: `blue.gradient`
```
0.0: 0.0, 0.0, 0.0
0.5: 0.0, 0.0, 1.0
1.0: 1.0, 1.0, 1.0
```

This defines a gradient that starts black (at position 0.0), is 100% blue in the middle (0.5) and turns into white at the end (1.0). The color at position 0.25 would be `R=0, G=0, B=0.5` and the color at 0.75 would be `R=0.5, G=0.5, B=1.0`.

The default colors at position 0.0 and 1.0 are black and white.

Points that are inside the Mandelbrot Set are always black.
