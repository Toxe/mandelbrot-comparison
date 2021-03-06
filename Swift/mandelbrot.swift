/*
 * Swift Mandelbrot.
 *
 * Usage: ./build/Swift/mandelbrot_swift <image width> <image height> <max iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
 * Example: ./build/Swift/mandelbrot_swift 800 600 200 1 -0.8 0.0 2.2 gradients/benchmark.gradient mandelbrot.raw
 *
 * Tobias Brückner, 2016
 */

import Foundation

enum RuntimeError: Error {
    case error(reason: String)
}

struct PixelColor {
    var r, g, b: UInt8
}

struct GradientColor {
    var pos: Float
    var r, g, b: Float
}

struct Gradient {
    var colors = [GradientColor]()
}

struct CalculationResult {
    var iter: Int
    var distanceToNextIteration: Float
}

func lerp(_ a: Float, _ b: Float, _ t: Float) -> Float {
    return (1.0 - t) * a + t * b
}

func lerp(_ a: Double, _ b: Double, _ t: Double) -> Double {
    return (1.0 - t) * a + t * b
}

func cumsum(_ iterationsHistogram: [Int32]) -> [Int32] {
    var cdf = [Int32](repeating: 0, count: iterationsHistogram.count)
    var total: Int32 = 0

    for (i, n) in iterationsHistogram.enumerated() {
        total += n
        cdf[i] = total
    }

    return cdf
}

// Compare two float values for "enough" equality.
func equalEnough(_ a: Float, _ b: Float) -> Bool {
    let absA = abs(a)
    let absB = abs(b)

    return abs(absA - absB) <= max(absA, absB) * Float.ulpOfOne
}

func gradientGetIndexOfColorAtPosition(_ gradient: Gradient, _ pos: Float) -> Int? {
    for i in 0..<gradient.colors.count {
        if equalEnough(gradient.colors[i].pos, pos) {
            return i
        }
    }

    return nil
}

func loadGradient(_ filename: String) throws -> Gradient {
    var gradient = Gradient()

    gradient.colors.append(GradientColor(pos: 0.0, r: 0.0, g: 0.0, b: 0.0))
    gradient.colors.append(GradientColor(pos: 1.0, r: 1.0, g: 1.0, b: 1.0))

    let contentsOfFile = try String(contentsOfFile: filename)
    let lines = contentsOfFile.components(separatedBy: "\n")
    let regex = try NSRegularExpression(pattern: "([0-9]*\\.?[0-9]+):\\s*([0-9]*\\.?[0-9]+),\\s*([0-9]*\\.?[0-9]+),\\s*([0-9]*\\.?[0-9]+)$")

    for line in lines {
        let matches = regex.matches(in: line, options: [], range: NSMakeRange(0, line.count))

        if matches.count == 1 {
            let pos = Float((line as NSString).substring(with: matches[0].range(at: 1)))
            let r = Float((line as NSString).substring(with: matches[0].range(at: 2)))
            let g = Float((line as NSString).substring(with: matches[0].range(at: 3)))
            let b = Float((line as NSString).substring(with: matches[0].range(at: 4)))

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

    gradient.colors.sort { $0.pos < $1.pos }
    return gradient
}

func colorFromGradientRange(_ left: GradientColor, _ right: GradientColor, _ pos: Float) -> PixelColor {
    let relativePosBetweenColors = (pos - left.pos) / (right.pos - left.pos)
    let r = lerp(left.r, right.r, relativePosBetweenColors)
    let g = lerp(left.g, right.g, relativePosBetweenColors)
    let b = lerp(left.b, right.b, relativePosBetweenColors)
    return PixelColor(r: UInt8(255.0 * r), g: UInt8(255.0 * g), b: UInt8(255.0 * b))
}

func colorFromGradient(_ gradient: Gradient, _ posInGradient: Float) -> PixelColor {
    var left = 0
    let colors = gradient.colors

    for right in 1 ..< colors.count {
        if posInGradient >= colors[left].pos && posInGradient <= colors[right].pos {
            return colorFromGradientRange(colors[left], colors[right], posInGradient)
        }

        left = right
    }

    return PixelColor(r: 0, g: 0, b: 0)
}

func mandelbrotCalc(_ imageWidth: Int, _ imageHeight: Int, _ maxIterations: Int, _ centerX: Double, _ centerY: Double, _ height: Double, _ iterationsHistogram: inout [Int32], _ resultsPerPoint: inout [CalculationResult]) {
    let width = height * (Double(imageWidth) / Double(imageHeight))

    let xLeft   = centerX - width / 2.0
    let xRight  = centerX + width / 2.0
    let yTop    = centerY + height / 2.0
    let yBottom = centerY - height / 2.0

    let bailout = 20.0
    let bailoutSquared = bailout * bailout
    let logLogBailout = log(log(bailout))
    let log2 = log(2.0)

    var finalMagnitude = 0.0

    for i in 0 ..< iterationsHistogram.count {
        iterationsHistogram[i] = 0
    }

    var pixel = 0

    for pixelY in 0 ..< imageHeight {
        let y0 = lerp(yTop, yBottom, Double(pixelY) / Double(imageHeight))

        for pixelX in 0 ..< imageWidth {
            let x0 = lerp(xLeft, xRight, Double(pixelX) / Double(imageWidth))

            var x = 0.0
            var y = 0.0

            // iteration, will be from 1 .. maxIterations once the loop is done
            var iter = 0

            while iter < maxIterations {
                let xSquared = x*x
                let ySquared = y*y

                if xSquared + ySquared >= bailoutSquared {
                    finalMagnitude = sqrt(xSquared + ySquared)
                    break
                }

                y = 2.0*x*y + y0
                x = xSquared - ySquared + x0

                iter += 1
            }

            if iter < maxIterations {
                iterationsHistogram[iter] += 1  // iter: 1 .. maxIterations-1, no need to count iterationsHistogram[max_iterations]
                resultsPerPoint[pixel] = CalculationResult(iter: iter, distanceToNextIteration: 1.0 - Float(min(1.0, (log(log(finalMagnitude)) - logLogBailout) / log2)))
            } else {
                resultsPerPoint[pixel] = CalculationResult(iter: 0, distanceToNextIteration: 0.0)
            }

            pixel += 1
        }
    }
}

func equalizeHistogram(_ iterationsHistogram: [Int32], _ maxIterations: Int) -> [Float] {
    // Calculate the CDF (Cumulative Distribution Function) by accumulating all iteration counts.
    // Element [0] is unused and iterationsHistogram[maxIterations] should be zero (as we do not count
    // the iterations of the points inside the Mandelbrot Set).
    let cdf = cumsum(iterationsHistogram)

    // Get the minimum value in the CDF that is bigger than zero and the sum of all iteration counts
    // from iterationsHistogram (which is the last value of the CDF).
    let cdfMin = cdf.first(where: { $0 > 0 })
    let totalIterations = cdf.last

    // normalize all values from the CDF that are bigger than zero to a range of 0.0 .. maxIterations
    let f = Float(maxIterations) / Float(totalIterations! - cdfMin!)
    return cdf.map { $0 > 0 ? f * Float($0 - cdfMin!) : 0.0 }
}

func mandelbrotColorize(_ imageWidth: Int, _ imageHeight: Int, _ maxIterations: Int, _ gradient: Gradient, _ imageData: inout [PixelColor], _ iterationsHistogram: [Int32], _ resultsPerPoint: [CalculationResult]) {
    let equalizedIterations = equalizeHistogram(iterationsHistogram, maxIterations)

    for (pixel, calcResult) in resultsPerPoint.enumerated() {
        if calcResult.iter == maxIterations {
            // points inside the Mandelbrot Set are always painted black
            imageData[pixel] = PixelColor(r: 0, g: 0, b: 0)
        } else {
            // The equalized iteration value (in the range of 0 .. maxIterations) represents the
            // position of the pixel color in the color gradiant and needs to be mapped to 0.0 .. 1.0.
            // To achieve smooth coloring we need to edge the equalized iteration towards the next
            // iteration, determined by the distance between the two iterations.
            let iterCurr = equalizedIterations[calcResult.iter]
            let iterNext = equalizedIterations[calcResult.iter + 1]

            let smoothedIteration = lerp(iterCurr, iterNext, calcResult.distanceToNextIteration)
            let posInGradient = smoothedIteration / Float(maxIterations)

            imageData[pixel] = colorFromGradient(gradient, posInGradient)
        }
    }
}

func saveImageData(_ filename: String, _ imageData: [PixelColor]) throws {
    let d = Data.init(bytes: imageData, count: MemoryLayout<PixelColor>.size * imageData.count)
    try d.write(to: URL(fileURLWithPath: filename, isDirectory: false))
}

func mean(_ values: [Double]) -> Double {
    return values.reduce(0.0, { $0 + $1 }) / Double(values.count)
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

func evalIntArgument(_ s: String, min: Int, max: Int) throws -> Int {
    guard let value = Int(s) else {
        throw RuntimeError.error(reason: "invalid value \(s)")
    }

    return value
}

func evalDoubleArgument(_ s: String, min: Double, max: Double) throws -> Double {
    guard let value = Double(s) else {
        throw RuntimeError.error(reason: "invalid value \(s)")
    }

    return value
}

func evalArguments(_ arguments: [String]) throws -> (Int, Int, Int, Int, Double, Double, Double, String, String)
{
    if arguments.count != 10 {
        throw RuntimeError.error(reason: "invalid number of arguments")
    }

    let imageWidth  = try evalIntArgument(arguments[1], min: 1, max: 100_000)
    let imageHeight = try evalIntArgument(arguments[2], min: 1, max: 100_000)
    let iterations  = try evalIntArgument(arguments[3], min: 1, max: 1_000_000_000)
    let repetitions = try evalIntArgument(arguments[4], min: 1, max: 1_000_000)
    let centerX     = try evalDoubleArgument(arguments[5], min: -100.0, max: 100.0)
    let centerY     = try evalDoubleArgument(arguments[6], min: -100.0, max: 100.0)
    let height      = try evalDoubleArgument(arguments[7], min: -100.0, max: 100.0)
    let colors      = arguments[8]
    let filename    = arguments[9]

    return (imageWidth, imageHeight, iterations, repetitions, centerX, centerY, height, colors, filename)
}

func go(_ imageWidth: Int, _ imageHeight: Int, _ maxIterations: Int, _ centerX: Double, _ centerY: Double, _ height: Double, _ gradient: Gradient, _ repetitions: Int) -> ([PixelColor], [Double]) {
    // iterationsHistogram: for simplicity we only use indices [1] .. [max_iterations], [0] is unused
    var iterationsHistogram = [Int32](repeating: 0, count: maxIterations + 1)

    // For every point store a tuple consisting of the final iteration and (for escaped points)
    // the distance to the next iteration (as value of 0.0 .. 1.0).
    var resultsPerPoint = [CalculationResult](repeating: CalculationResult(iter: 0, distanceToNextIteration: 0.0), count: imageWidth * imageHeight)

    var imageData = [PixelColor](repeating: PixelColor(r: 0, g: 0, b: 0), count: imageWidth * imageHeight)
    var durations = [Double]()

    for _ in 0 ..< repetitions {
        let t1 = Date()
        mandelbrotCalc(imageWidth, imageHeight, maxIterations, centerX, centerY, height, &iterationsHistogram, &resultsPerPoint)
        mandelbrotColorize(imageWidth, imageHeight, maxIterations, gradient, &imageData, iterationsHistogram, resultsPerPoint)
        let t2 = Date()

        durations.append(t2.timeIntervalSince(t1))
    }

    return (imageData, durations)
}

func main() throws {
    let (imageWidth, imageHeight, iterations, repetitions, centerX, centerY, height, gradientFilename, filename) = try evalArguments(CommandLine.arguments)
    let gradient = try loadGradient(gradientFilename)

    let (imageData, durations) = go(imageWidth, imageHeight, iterations, centerX, centerY, height, gradient, repetitions)

    try saveImageData(filename, imageData)
    showSummary(durations)
}
