/*
 * Swift Mandelbrot.
 *
 * Usage: mandelbrot_swift <image_width> <image_height> <iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
 * Example: mandelbrot_swift 320 200 20 1 -0.5 0.0 2.0 blue.gradient mandelbrot.raw
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

struct GradientColor {
    var pos: Double
    var r, g, b: Double
}

struct Gradient {
    var colors = [GradientColor]()
}


// Compare two double values for "enough" equality.
func equalEnough(_ a: Double, _ b: Double) -> Bool {
    let absA = abs(a)
    let absB = abs(b)

    return abs(absA - absB) <= max(absA, absB) * Double.ulpOfOne
}

func gradientGetIndexOfColorAtPosition(_ gradient: Gradient, _ pos: Double) -> Int? {
    for i in 0..<gradient.colors.count {
        if equalEnough(gradient.colors[i].pos, pos) {
            return i
        }
    }

    return nil
}

func loadGradient(_ filename: String) -> Gradient? {
    var gradient = Gradient()

    gradient.colors.append(GradientColor(pos: 0.0, r: 0.0, g: 0.0, b: 0.0))
    gradient.colors.append(GradientColor(pos: 1.0, r: 1.0, g: 1.0, b: 1.0))

    do {
        let contentsOfFile = try String(contentsOfFile: filename)
        let lines = contentsOfFile.components(separatedBy: "\n")
        let regex = try NSRegularExpression(pattern: "([0-9]*\\.?[0-9]+):\\s*([0-9]*\\.?[0-9]+),\\s*([0-9]*\\.?[0-9]+),\\s*([0-9]*\\.?[0-9]+)")

        for line in lines {
            let matches = regex.matches(in: line, options: [], range: NSMakeRange(0, line.count))

            if matches.count == 1 {
                if matches[0].numberOfRanges == 4+1 {
                    let pos = Double((line as NSString).substring(with: matches[0].range(at: 1)))
                    let r = Double((line as NSString).substring(with: matches[0].range(at: 2)))
                    let g = Double((line as NSString).substring(with: matches[0].range(at: 3)))
                    let b = Double((line as NSString).substring(with: matches[0].range(at: 4)))

                    if pos != nil && r != nil && g != nil && b != nil {
                        if let index = gradientGetIndexOfColorAtPosition(gradient, pos!) {
                            gradient.colors[index].r = r!
                            gradient.colors[index].g = g!
                            gradient.colors[index].b = b!
                        } else {
                            gradient.colors.append(GradientColor(pos: pos!, r: r!, g: g!, b: b!))
                        }
                    }
                }
            }
        }
    } catch {
        return nil
    }

    gradient.colors.sort { $0.pos < $1.pos }

    return gradient
}

func colorFromGradientRange(_ left: GradientColor, _ right: GradientColor, _ pos: Double) -> (Double, Double, Double) {
    let pos2 = (pos - left.pos) / (right.pos - left.pos)

    return (((right.r - left.r) * pos2) + left.r,
            ((right.g - left.g) * pos2) + left.g,
            ((right.b - left.b) * pos2) + left.b)
}

func colorFromGradient(_ gradient: Gradient, _ posInGradient: Double) -> (Double, Double, Double) {
    var left = 0
    let colors = gradient.colors

    for right in 1 ..< colors.count {
        if posInGradient >= colors[left].pos && posInGradient <= colors[right].pos {
            return colorFromGradientRange(colors[left], colors[right], posInGradient)
        }

        left = right
    }

    return (0.0, 0.0, 0.0)
}

func mandelbrotCalc(_ imageWidth: Int, _ imageHeight: Int, _ maxIterations: Int, _ centerX: Double, _ centerY: Double, _ height: Double, _ histogram: inout [Int], _ iterationsPerPixel: inout [Int], _ smoothedDistancesToNextIterationPerPixel: inout [Double]) {
    let width = height * (Double(imageWidth) / Double(imageHeight))

    let xLeft   = centerX - width / 2.0
//  let xRight  = centerX + width / 2.0
    let yTop    = centerY + height / 2.0
//  let yBottom = centerY - height / 2.0

    let bailout = 20.0
    let bailoutSquared = bailout * bailout
    let logLogBailout = log(log(bailout))
    let log2 = log(2.0)

    var finalMagnitude = 0.0

    for i in 0 ..< histogram.count {
        histogram[i] = 0
    }

    for pixelY in 0 ..< imageHeight {
        let y0 = yTop - height * (Double(pixelY) / Double(imageHeight))

        for pixelX in 0 ..< imageWidth {
            let x0 = xLeft + width * (Double(pixelX) / Double(imageWidth))

            var x = 0.0
            var y = 0.0

            // iteration, will be from 1 to max_iterations once the loop is done
            var iter = 0

            while iter < maxIterations {
                let xSquared = x*x
                let ySquared = y*y

                if xSquared + ySquared >= bailoutSquared {
                    finalMagnitude = sqrt(xSquared + ySquared)
                    break
                }

                let xtemp = xSquared - ySquared + x0
                y = 2.0*x*y + y0
                x = xtemp

                iter += 1
            }

            let pixel = pixelY * imageWidth + pixelX

            if iter < maxIterations {
                smoothedDistancesToNextIterationPerPixel[pixel] = 1.0 - min(1.0, (log(log(finalMagnitude)) - logLogBailout) / log2)
                histogram[iter] += 1  // no need to count histogram[max_iterations]
            }

            iterationsPerPixel[pixel] = iter  // 1 .. max_iterations
        }
    }
}

func mandelbrotColorize(_ imageWidth: Int, _ imageHeight: Int, _ maxIterations: Int, _ gradient: Gradient, _ imageData: inout [UInt8], _ histogram: [Int], _ iterationsPerPixel: [Int], _ smoothedDistancesToNextIterationPerPixel: [Double], _ normalizedColors: inout [Double]) {
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
            let pixel = pixelY * imageWidth + pixelX
            let iter = iterationsPerPixel[pixel]  // 1 .. max_iterations

            if iter == maxIterations {
                // pixels with max. iterations (aka. inside the Mandelbrot Set) are always black
                imageData[3 * pixel + 0] = 0
                imageData[3 * pixel + 1] = 0
                imageData[3 * pixel + 2] = 0
            } else {
                // we use the color of the previous iteration in order to cover the full gradient range
                let colorOfPreviousIter = normalizedColors[iter - 1]
                let colorOfCurrentIter  = normalizedColors[iter]

                let smoothedDistanceToNextIteration = smoothedDistancesToNextIterationPerPixel[pixel]  // 0 .. <1.0
                let posInGradient = colorOfPreviousIter + smoothedDistanceToNextIteration * (colorOfCurrentIter - colorOfPreviousIter)

                let (r, g, b) = colorFromGradient(gradient, posInGradient)

                imageData[3 * pixel + 0] = UInt8(255.0 * r)
                imageData[3 * pixel + 1] = UInt8(255.0 * g)
                imageData[3 * pixel + 2] = UInt8(255.0 * b)
            }
        }
    }
}

func saveImageData(_ filename: String, _ imageData: [UInt8]) -> Bool {
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

    for value in values {
        sum += value
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

func evalIntArgument(_ s: String, min: Int, max: Int) -> Int {
    guard let value = Int(s) else {
        die(ExitCode.EVAL_ARGS)
    }

    return value
}

func evalDoubleArgument(_ s: String, min: Double, max: Double) -> Double {
    guard let value = Double(s) else {
        die(ExitCode.EVAL_ARGS)
    }

    return value
}

func evalArguments(_ arguments: [String]) -> (Int, Int, Int, Int, Double, Double, Double, String, String)
{
    if arguments.count < 10 {
        die(ExitCode.EVAL_ARGS)
    }

    let imageWidth  = evalIntArgument(arguments[1], min: 1, max: 100_000)
    let imageHeight = evalIntArgument(arguments[2], min: 1, max: 100_000)
    let iterations  = evalIntArgument(arguments[3], min: 1, max: 1_000_000_000)
    let repetitions = evalIntArgument(arguments[4], min: 1, max: 1_000_000)
    let centerX     = evalDoubleArgument(arguments[5], min: -100.0, max: 100.0)
    let centerY     = evalDoubleArgument(arguments[6], min: -100.0, max: 100.0)
    let height      = evalDoubleArgument(arguments[7], min: -100.0, max: 100.0)
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
    let (imageWidth, imageHeight, iterations, repetitions, centerX, centerY, height, gradientFilename, filename) = evalArguments(CommandLine.arguments)

    guard let gradient = loadGradient(gradientFilename) else {
        die(ExitCode.LOAD_GRADIENT)
    }

    var imageData = [UInt8](repeating: 0, count: imageWidth * imageHeight * 3)
    var durations = [Double]()

    go(imageWidth, imageHeight, iterations, centerX, centerY, height, gradient, &imageData, &durations, repetitions)

    if !saveImageData(filename, imageData) {
        die(ExitCode.SAVE_IMAGE)
    }

    showSummary(durations)
}
