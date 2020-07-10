version development

workflow SepWorkflow2 {
  input {}
  call SepTest
}

task SepTest {
    input {
        Array[String] inp = ["value1", "value2", "value3"]
    }
    command {}
    output {
        String out = sep(",", inp)
    }
}