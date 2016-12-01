workflow square {
    Array[Int] arr = [1, 2]
    
    scatter(i in arr) {
        Int double = i * i
    }
    
    output {
        Array[Int] doubled = double
    }
}
