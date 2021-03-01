do {
    try main()
} catch RuntimeError.error(let reason) {
    print("Runtime error: \(reason)")
} catch {
    print("Runtime error: \(error)")
}
