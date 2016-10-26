/*
 * Swift Mandelbrot.
 *
 * Usage: mandelbrot <image_width> <image_height> <iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
 * Example: mandelbrot 320 200 20 1 -0.5 0.0 2.0 blue.gradient mandelbrot.raw
 *
 * Gradient file example:
 *   0.0: 0.0, 0.0, 0.0
 *   0.5: 0.0, 0.0, 1.0
 *   1.0: 1.0, 1.0, 1.0
 *
 * Tobias Brückner, 2016
 */

import Foundation

enum ExitCode: Int32 {
    case ALLOC_MEMORY = 1
    case EVAL_ARGS
    case LOAD_GRADIENT
    case SAVE_IMAGE
    case GETTIME
}

class GradientColor {
    var pos: Double
    var r, g, b: Double

    init(_ pos: Double, _ r: Double, _ g: Double, _ b: Double) {
        self.pos = pos
        self.r = r
        self.g = g
        self.b = b
    }
}

class Gradient {
    var colors = [GradientColor]()
}


// Compare two double values for "enough" equality.
func equalEnough(_ a: Double, _ b: Double) -> Bool {
    let absA = abs(a)
    let absB = abs(b)

    return abs(absA - absB) <= max(absA, absB) * DBL_EPSILON
}

func gradientGetColorAtPosition(_ gradient: Gradient, _ pos: Double) -> GradientColor? {
    for i in 0 ..< gradient.colors.count {
        if equalEnough(gradient.colors[i].pos, pos) {
            return gradient.colors[i]
        }
    }

    return nil
}

func cmpColorPosFunc(_ a: GradientColor, _ b: GradientColor) -> Bool {
    return a.pos < b.pos
}

func loadGradient(_ filename: String) -> Gradient? {
    let gradient = Gradient()

    gradient.colors.append(GradientColor(0.0, 0.0, 0.0, 0.0))
    gradient.colors.append(GradientColor(1.0, 1.0, 1.0, 1.0))

    do {
        let contentsOfFile = try String(contentsOfFile: filename)
        let lines = contentsOfFile.components(separatedBy: "\n")
        let regex = try NSRegularExpression(pattern: "([0-9]*\\.?[0-9]+):\\s*([0-9]*\\.?[0-9]+),\\s*([0-9]*\\.?[0-9]+),\\s*([0-9]*\\.?[0-9]+)")

        for i in 0 ..< lines.count {
            let line = lines[i]
            let matches = regex.matches(in: line, options: [], range: NSMakeRange(0, line.characters.count))

            if matches.count == 1 {
                if matches[0].numberOfRanges == 4+1 {
                    let pos = Double((line as NSString).substring(with: matches[0].rangeAt(1)))
                    let r = Double((line as NSString).substring(with: matches[0].rangeAt(2)))
                    let g = Double((line as NSString).substring(with: matches[0].rangeAt(3)))
                    let b = Double((line as NSString).substring(with: matches[0].rangeAt(4)))

                    if pos != nil && r != nil && g != nil && b != nil {
                        let color = gradientGetColorAtPosition(gradient, pos!)

                        if color != nil {
                            color!.r = r!
                            color!.g = g!
                            color!.b = b!
                        } else {
                            gradient.colors.append(GradientColor(pos!, r!, g!, b!))
                        }
                    }
                }
            }
        }
    } catch {
        return nil
    }

    gradient.colors.sort(by: cmpColorPosFunc)
    return gradient
}

func colorFromGradientRange(_ left: GradientColor, _ right: GradientColor, _ pos: Double, _ r: inout Double, _ g: inout Double, _ b: inout Double) {
    let pos2 = (pos - left.pos) / (right.pos - left.pos)

    r = ((right.r - left.r) * pos2) + left.r
    g = ((right.g - left.g) * pos2) + left.g
    b = ((right.b - left.b) * pos2) + left.b
}

func colorFromGradient(_ gradient: Gradient, _ posInGradient: Double, _ r: inout Double, _ g: inout Double, _ b: inout Double) -> Bool {
    var pos = posInGradient

    if (pos < 0.0) {
        pos = 0.0
    }

    if (pos > 1.0) {
        pos = 1.0
    }

    var left, right: GradientColor
    left = gradient.colors[0]

    for i in 1 ..< gradient.colors.count {
        right = gradient.colors[i]

        if pos >= left.pos && pos <= right.pos {
            colorFromGradientRange(left, right, pos, &r, &g, &b)
            return true
        }

        left = right
    }

    return false
}

func mandelbrotCalc(_ imageWidth: Int, _ imageHeight: Int, _ maxIterations: Int, _ centerX: Double, _ centerY: Double, _ height: Double, _ histogram: inout [Int], _ iterationsPerPixel: inout [Int], _ smoothedDistancesToNextIterationPerPixel: inout [Double]) {
    let width = height * (Double(imageWidth) / Double(imageHeight))

    let xLeft   = centerX - width / 2.0
//  let xRight  = centerX + width / 2.0
    let yTop    = centerY + height / 2.0
//  let yBottom = centerY - height / 2.0

    var iter: Int
    var x0, y0: Double
    var x, y: Double
    var xtemp: Double
    var xSquared = 0.0, ySquared = 0.0

    let bailout = 20.0
    let bailoutSquared = bailout * bailout
    let logLogBailout = log(log(bailout))
    let log2 = log(2.0)

    for i in 0 ..< histogram.count {
        histogram[i] = 0
    }

    for pixelY in 0 ..< imageHeight {
        for pixelX in 0 ..< imageWidth {
            x0 = xLeft + width * (Double(pixelX) / Double(imageWidth))
            y0 = yTop - height * (Double(pixelY) / Double(imageHeight))

            x = 0.0
            y = 0.0

            // iteration, will be from 1 to max_iterations once the loop is done
            iter = 0

            while iter < maxIterations {
                xSquared = x*x
                ySquared = y*y

                if xSquared + ySquared >= bailoutSquared {
                    break
                }

                xtemp = xSquared - ySquared + x0
                y = 2.0*x*y + y0
                x = xtemp

                iter += 1
            }

            if iter < maxIterations {
                let finalMagnitude = sqrt(xSquared + ySquared)
                smoothedDistancesToNextIterationPerPixel[pixelY * imageWidth + pixelX] = 1.0 - min(1.0, (log(log(finalMagnitude)) - logLogBailout) / log2)
            }

            histogram[iter] += 1
            iterationsPerPixel[pixelY * imageWidth + pixelX] = iter  // 1 .. max_iterations

        }
    }
}

func mandelbrotColorize(_ imageWidth: Int, _ imageHeight: Int, _ maxIterations: Int, _ gradient: Gradient, _ imageData: inout [UInt8], _ histogram: [Int], _ iterationsPerPixel: [Int], _ smoothedDistancesToNextIterationPerPixel: [Double], _ normalizedColors: inout [Double]) {
    for i in 0 ..< normalizedColors.count {
        normalizedColors[i] = 0.0
    }

    // Sum all iterations, not counting the last one at position histogram[max_iterations] (which
    // are points in the Mandelbrot Set).
    var totalIterations = 0

    for i in 1 ..< maxIterations {
        totalIterations += histogram[i]
    }

    // Normalize the colors (0.0 .. 1.0) based on how often they are used in the image, not counting
    // histogram[max_iterations] (which are points in the Mandelbrot Set).
    var runningTotal = 0

    for i in 1 ..< maxIterations {
        runningTotal += histogram[i]
        normalizedColors[i] = Double(runningTotal) / Double(totalIterations)
    }

    for pixelY in 0 ..< imageHeight {
        for pixelX in 0 ..< imageWidth {
            let iter = iterationsPerPixel[pixelY * imageWidth + pixelX]  // 1 .. max_iterations
            var r = 0.0, g = 0.0, b = 0.0

            if iter == maxIterations {
                // pixels with max. iterations (aka. inside the Mandelbrot Set) are always black
                r = 0.0
                g = 0.0
                b = 0.0
            } else {
                // we use the color of the previous iteration in order to cover the full gradient range
                let colorOfPreviousIter = normalizedColors[iter - 1]
                let colorOfCurrentIter  = normalizedColors[iter]

                let smoothedDistanceToNextIteration = smoothedDistancesToNextIterationPerPixel[pixelY * imageWidth + pixelX]  // 0 .. <1.0
                let posInGradient = colorOfPreviousIter + smoothedDistanceToNextIteration * (colorOfCurrentIter - colorOfPreviousIter)

                let _ = colorFromGradient(gradient, posInGradient, &r, &g, &b)
            }

            imageData[3 * (pixelY * imageWidth + pixelX) + 0] = UInt8(255.0 * r)
            imageData[3 * (pixelY * imageWidth + pixelX) + 1] = UInt8(255.0 * g)
            imageData[3 * (pixelY * imageWidth + pixelX) + 2] = UInt8(255.0 * b)
        }
    }
}

func saveImage(_ filename: String, _ imageData: [UInt8]) -> Bool {
    do {
        let d = Data(imageData)
        try d.write(to: URL(fileURLWithPath: filename, isDirectory: false))
    } catch {
        return false
    }

    return true
}

func mean(_ values: [Double]) -> Double {
    var sum = 0.0

    for i in 0 ..< values.count {
        sum += values[i]
    }

    return sum / Double(values.count)
}

func median(_ values: [Double]) -> Double {
    let sortedValues = values.sorted()

    if sortedValues.count % 2 == 1 {
        return sortedValues[(sortedValues.count - 1) / 2]
    } else {
        return (sortedValues[sortedValues.count/2 - 1] + sortedValues[sortedValues.count/2]) / 2.0
    }
}

func showSummary(_ durations: [Double]) {
    if durations.count == 1 {
        print("\(durations[0]) s")
    } else {
        print("mean: \(mean(durations)) s, median: \(median(durations)) s (repetitions=\(durations.count)) \(durations.sorted())")
    }
}

func die(_ error: ExitCode) -> Never {
    print("Error: \(error.rawValue)")
    exit(error.rawValue)
}

func evalIntArg(_ s: String, min: Int, max: Int) -> Int {
    let value = Int(s)

    if value == nil {
        die(ExitCode.EVAL_ARGS)
    }

    return value!
}

func evalDoubleArg(_ s: String, min: Double, max: Double) -> Double {
    let value = Double(s)

    if value == nil {
        die(ExitCode.EVAL_ARGS)
    }

    return value!
}

func evalArgs(_ arguments: [String]) -> (Int, Int, Int, Int, Double, Double, Double, String, String)
{
    if arguments.count < 10 {
        die(ExitCode.EVAL_ARGS)
    }

    let imageWidth  = evalIntArg(arguments[1], min: 1, max: 100_000)
    let imageHeight = evalIntArg(arguments[2], min: 1, max: 100_000)
    let iterations  = evalIntArg(arguments[3], min: 1, max: 1_000_000_000)
    let repetitions = evalIntArg(arguments[4], min: 1, max: 1_000_000)
    let centerX     = evalDoubleArg(arguments[5], min: -100.0, max: 100.0)
    let centerY     = evalDoubleArg(arguments[6], min: -100.0, max: 100.0)
    let height      = evalDoubleArg(arguments[7], min: -100.0, max: 100.0)
    let colors      = arguments[8]
    let filename    = arguments[9]

    return (imageWidth, imageHeight, iterations, repetitions, centerX, centerY, height, colors, filename)
}

func go(_ imageWidth: Int, _ imageHeight: Int, _ maxIterations: Int, _ centerX: Double, _ centerY: Double, _ height: Double, _ gradient: Gradient, _ imageData: inout [UInt8], _ durations: inout [Double], _ repetitions: Int) {
    // histogram & normalized_colors: for simplicity we only use indices [1] .. [max_iterations], [0] is unused
    var histogram = [Int](repeating: 0, count: maxIterations + 1)
    var normalizedColors = [Double](repeating: 0, count: maxIterations + 1)
    var iterationsPerPixel = [Int](repeating: 0, count: imageWidth * imageHeight)
    var smoothedDistancesToNextIterationPerPixel = [Double](repeating: 0, count: imageWidth * imageHeight)

    for _ in 0 ..< repetitions {
        let t1 = Date()
        mandelbrotCalc(imageWidth, imageHeight, maxIterations, centerX, centerY, height, &histogram, &iterationsPerPixel, &smoothedDistancesToNextIterationPerPixel)
        mandelbrotColorize(imageWidth, imageHeight, maxIterations, gradient, &imageData, histogram, iterationsPerPixel, smoothedDistancesToNextIterationPerPixel, &normalizedColors)
        let t2 = Date()

        durations.append(t2.timeIntervalSince(t1))
    }
}

func main() {
    let (imageWidth, imageHeight, iterations, repetitions, centerX, centerY, height, gradientFilename, filename) = evalArgs(CommandLine.arguments)

    let gradient = loadGradient(gradientFilename)

    if gradient == nil {
        die(ExitCode.LOAD_GRADIENT)
    }

    var imageData = [UInt8](repeating: 0, count: imageWidth * imageHeight * 3)
    var durations = [Double]()

    go(imageWidth, imageHeight, iterations, centerX, centerY, height, gradient!, &imageData, &durations, repetitions)

    if !saveImage(filename, imageData) {
        die(ExitCode.SAVE_IMAGE)
    }

    showSummary(durations)
}
