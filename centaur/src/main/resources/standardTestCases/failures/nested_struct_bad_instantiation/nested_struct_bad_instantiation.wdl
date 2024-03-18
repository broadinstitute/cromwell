version 1.0


struct firstLayer {
    String first
    Int number
}

struct secondLayer {
    String second
    firstLayer lowerLayer
}


workflow MinimalStructExample {
    input {
        String test
        Int integer
    }

    firstLayer example1 = {"first": test, "number": integer}
    secondLayer example2 = {"second": test, "lowerLayer": example1}

    output {
        String example3 = example2.second
    }
}
