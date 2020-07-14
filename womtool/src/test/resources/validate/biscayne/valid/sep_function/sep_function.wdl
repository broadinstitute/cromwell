version development

workflow SepWorkflow {
  input {}
  call SepTestInInterpolatorBlock
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