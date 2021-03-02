<?php
// PHP Mandelbrot.
//
// Usage: php -d xdebug.mode=off PHP/mandelbrot.php <image width> <image height> <max iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
// Example: php -d xdebug.mode=off PHP/mandelbrot.php 800 600 200 1 -0.8 0.0 2.2 gradients/benchmark.gradient mandelbrot.raw
//
// Tobias Brückner, 2021

class GradientColor
{
    public $pos;
    public $r, $g, $b;

    function __construct($pos, $r, $g, $b)
    {
        $this->pos = (float) $pos;
        $this->update_colors($r, $g, $b);
    }

    function update_colors($r, $g, $b)
    {
        $this->r = (float) $r;
        $this->g = (float) $g;
        $this->b = (float) $b;
    }
}

class Gradient
{
    public $colors = [];
}

function lerp($a, $b, $t)
{
    return (1.0 - $t) * $a + $t * $b;
}

function cumsum($hist)
{
    $cdf = [];
    $total = 0;

    foreach ($hist as $n) {
    	$total += $n;
    	$cdf[] = $total;
    }

    return $cdf;
}

// Compare two float values for "enough" equality.
function equal_enough($a, $b)
{
    $a = abs($a);
    $b = abs($b);
    return abs($a - $b) <= max($a, $b) * PHP_FLOAT_EPSILON;
}

function gradient_get_color_at_position($gradient, $pos)
{
    return current(array_filter($gradient->colors, fn($c) => equal_enough($c->pos, $pos)));
}

function load_gradient($gradient_filename)
{
    $re = '/([0-9]*\.?[0-9]+):\s*([0-9]*\.?[0-9]+),\s*([0-9]*\.?[0-9]+),\s*([0-9]*\.?[0-9]+)$/';

    $gradient = new Gradient();
    $gradient->colors[] = new GradientColor(0.0, 0.0, 0.0, 0.0);
    $gradient->colors[] = new GradientColor(1.0, 1.0, 1.0, 1.0);

    if (!($fp = fopen($gradient_filename, "r")))
        throw new RuntimeException("unable to open gradient file");

    while ($line = fgets($fp)) {
        if (preg_match($re, $line, $m)) {
            $col = gradient_get_color_at_position($gradient, floatval($m[1]));

            if ($col)
                $col->update_colors($m[2], $m[3], $m[4]);
            else
                $gradient->colors[] = new GradientColor($m[1], $m[2], $m[3], $m[4]);
        }
    }

    fclose($fp);

    usort($gradient->colors, fn($a, $b) => $a->pos <=> $b->pos);
    return $gradient;
}

function color_from_gradient_range($left_color, $right_color, $pos)
{
    $relative_pos_between_colors = ($pos - $left_color->pos) / ($right_color->pos - $left_color->pos);
    $r = lerp($left_color->r, $right_color->r, $relative_pos_between_colors);
    $g = lerp($left_color->g, $right_color->g, $relative_pos_between_colors);
    $b = lerp($left_color->b, $right_color->b, $relative_pos_between_colors);
    return [$r, $g, $b];
}

function color_from_gradient($gradient, $pos)
{
    $pairs = array_map(NULL, array_slice($gradient->colors, 0, -1), array_slice($gradient->colors, 1));

    foreach ($pairs as $colors)
        if ($pos >= $colors[0]->pos && $pos <= $colors[1]->pos)
            return color_from_gradient_range($colors[0], $colors[1], $pos);

    return NULL;
}

function mandelbrot_calc($image_width, $image_height, $max_iterations, $center_x, $center_y, $height, &$iterations_per_pixel, &$distances_to_next_iteration_per_pixel)
{
    $width = $height * ($image_width / $image_height);

    $x_left   = $center_x - $width / 2.0;
    $x_right  = $center_x + $width / 2.0;
    $y_top    = $center_y + $height / 2.0;
    $y_bottom = $center_y - $height / 2.0;

    $bailout = 20.0;
    $bailout_squared = $bailout * $bailout;
    $log_log_bailout = log(log($bailout));
    $log_2 = log(2.0);

    // for simplicity we only use indices [1] .. [max_iterations], [0] is unused
    $iterations_histogram = array_fill(0, $max_iterations + 1, 0);

    $pixel = 0;

    for ($pixel_y = 0; $pixel_y < $image_height; ++$pixel_y) {
        $y0 = lerp($y_top, $y_bottom, $pixel_y / $image_height);

        for ($pixel_x = 0; $pixel_x < $image_width; ++$pixel_x) {
            $x0 = lerp($x_left, $x_right, $pixel_x / $image_width);

            $x = 0.0;
            $y = 0.0;

            $x_squared = 0.0;
            $y_squared = 0.0;

            // iteration, will be from 1 .. max_iterations once the loop is done
            $iter = 0;

            while ($iter < $max_iterations) {
                $x_squared = $x * $x;
                $y_squared = $y * $y;

                if ($x_squared + $y_squared >= $bailout_squared)
                    break;

                $y = 2.0 * $x * $y + $y0;
                $x = $x_squared - $y_squared + $x0;

                ++$iter;
            }

            if ($iter < $max_iterations) {
                $final_magnitude = sqrt($x_squared + $y_squared);
                $distances_to_next_iteration_per_pixel[$pixel] = 1.0 - min(1.0, (log(log($final_magnitude)) - $log_log_bailout) / $log_2);
                $iterations_histogram[$iter] += 1;  // iter: 1 .. max_iterations-1, no need to count iterations_histogram[max_iterations]
            }

            $iterations_per_pixel[$pixel] = $iter;
            ++$pixel;
        }
    }

    return $iterations_histogram;
}

function equalize_histogram($iterations_histogram, $max_iterations)
{
    // Calculate the CDF (Cumulative Distribution Function) by accumulating all iteration counts.
    // Element [0] is unused and iterations_histogram[max_iterations] should be zero (as we do not count
    // the iterations of the points inside the Mandelbrot Set).
    $cdf = cumsum($iterations_histogram);

    // Get the minimum value in the CDF that is bigger than zero and the sum of all iteration counts
    // from iterations_histogram (which is the last value of the CDF).
    $cdf_min = current(array_filter($cdf, fn($x) => $x > 0));
    $total_iterations = $cdf[count($cdf) - 1];

    // normalize all values from the CDF that are bigger than zero to a range of 0.0 .. max_iterations
    $f = $max_iterations / ($total_iterations - $cdf_min);
    return array_map(fn($c) => $c > 0 ? $f * ($c - $cdf_min) : 0, $cdf);
}

function mandelbrot_colorize($image_width, $image_height, $max_iterations, $gradient, &$image_data, $iterations_histogram, $iterations_per_pixel, $distances_to_next_iteration_per_pixel)
{
    $equalized_iterations = equalize_histogram($iterations_histogram, $max_iterations);

    for ($pixel = 0; $pixel < $image_width * $image_height; ++$pixel) {
        $iter = $iterations_per_pixel[$pixel];  // 1 .. max_iterations

        if ($iter == $max_iterations) {
            // points inside the Mandelbrot Set are always painted black
            $image_data[3 * $pixel + 0] = 0;
            $image_data[3 * $pixel + 1] = 0;
            $image_data[3 * $pixel + 2] = 0;
        } else {
            // The equalized iteration value (in the range of 0 .. max_iterations) represents the
            // position of the pixel color in the color gradiant and needs to be mapped to 0.0 .. 1.0.
            // To achieve smooth coloring we need to edge the equalized iteration towards the next
            // iteration, determined by the distance between the two iterations.
            $iter_curr = $equalized_iterations[$iter];
            $iter_next = $equalized_iterations[$iter + 1];

            $smoothed_iteration = lerp($iter_curr, $iter_next, $distances_to_next_iteration_per_pixel[$pixel]);
            $pos_in_gradient = $smoothed_iteration / $max_iterations;

            [$r, $g, $b] = color_from_gradient($gradient, $pos_in_gradient);

            $image_data[3 * $pixel + 0] = intval(255.0 * $r);
            $image_data[3 * $pixel + 1] = intval(255.0 * $g);
            $image_data[3 * $pixel + 2] = intval(255.0 * $b);
        }
    }
}

function save_image($filename, $image_data)
{
    if (!($fp = fopen($filename, "wb")))
        throw new RuntimeException("unable to open output file");

    fwrite($fp, pack("C*", ...$image_data));
    fclose($fp);
}

function mean($values)
{
    return array_sum($values) / count($values);
}

function median($values)
{
    $sorted_values = $values;
    sort($sorted_values);

    if (count($sorted_values) % 2)
        return $sorted_values[(count($sorted_values) - 1) / 2];
    else
        return ($sorted_values[count($sorted_values) / 2 - 1] + $sorted_values[count($sorted_values) / 2]) / 2.0;
}

function show_summary($durations)
{
    if (count($durations) == 1) {
        printf("%f s\n", $durations[0]);
    } else {
        $sorted_values = $durations;
        sort($sorted_values);
        printf("mean: %f s, median: %f s (repetitions=%d) [%s]\n", mean($sorted_values), median($sorted_values), count($sorted_values), implode(", ", $sorted_values));
    }
}

function eval_int_arg($s, $min, $max)
{
    $value = intval($s);

    if ($value < $min || $value > $max)
        throw new RuntimeException("invalid value {$s}");

    return $value;
}

function eval_float_arg($s, $min, $max)
{
    $value = floatval($s);

    if ($value < $min || $value > $max || is_nan($value) || is_infinite(($value)))
        throw new RuntimeException("invalid value {$s}");

    return $value;
}

function eval_args()
{
    if ($_SERVER["argc"] != 10)
        throw new RuntimeException("invalid number of arguments");

    $image_width       = eval_int_arg($_SERVER["argv"][1], 1, 100000);
    $image_height      = eval_int_arg($_SERVER["argv"][2], 1, 100000);
    $max_iterations    = eval_int_arg($_SERVER["argv"][3], 1, 1000000000);
    $repetitions       = eval_int_arg($_SERVER["argv"][4], 1, 1000000);
    $center_x          = eval_float_arg($_SERVER["argv"][5], -100.0, 100.0);
    $center_y          = eval_float_arg($_SERVER["argv"][6], -100.0, 100.0);
    $height            = eval_float_arg($_SERVER["argv"][7], -100.0, 100.0);
    $gradient_filename = $_SERVER["argv"][8];
    $filename          = $_SERVER["argv"][9];

    return [$image_width, $image_height, $max_iterations, $repetitions, $center_x, $center_y, $height, $gradient_filename, $filename];
}

function go($image_width, $image_height, $max_iterations, $center_x, $center_y, $height, $gradient, &$image_data, $repetitions)
{
    $iterations_per_pixel = array_fill(0, $image_width * $image_height, 0);
    $distances_to_next_iteration_per_pixel = array_fill(0, $image_width * $image_height, 0);
    $durations = [];

    for ($i = 0; $i < $repetitions; ++$i) {
        $t1 = microtime(true);
        $iterations_histogram = mandelbrot_calc($image_width, $image_height, $max_iterations, $center_x, $center_y, $height, $iterations_per_pixel, $distances_to_next_iteration_per_pixel);
        mandelbrot_colorize($image_width, $image_height, $max_iterations, $gradient, $image_data, $iterations_histogram, $iterations_per_pixel, $distances_to_next_iteration_per_pixel);
        $t2 = microtime(true);

        $durations[] = $t2 - $t1;
    }

    return $durations;
}

function main()
{
    [$image_width, $image_height, $max_iterations, $repetitions, $center_x, $center_y, $height, $gradient_filename, $filename] = eval_args();
    $gradient = load_gradient($gradient_filename);
    $image_data = array_fill(0, 3 * $image_width * $image_height, 0);

    $durations = go($image_width, $image_height, $max_iterations, $center_x, $center_y, $height, $gradient, $image_data, $repetitions);

    save_image($filename, $image_data);
    show_summary($durations);
}

main();
