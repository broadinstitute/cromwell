version 1.0

task SepTest {
    input {
        Array[String] inp = ["value1", "value2", "value3"]
    }
    command {}
    output {
        String out = sep(",", inp)
    }
}

task SepTestInInterpolatorBlock {
    input {
        Array[String] inp = ["value1", "value2", "value3"]
    }
    command <<<
        echo ~{sep(",", inp)}
    >>>

    output {
        String out = read_string(stdout())
    }
}