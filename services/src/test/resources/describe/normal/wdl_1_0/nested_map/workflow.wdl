version 1.0

workflow wdl_1_0_map {

  input {
    Map[Int, Map[Map[String, Int], String]] input_1
  }

  output {
     Map[Int, Map[Map[String, Int], String]] output_1 = input_1
  }

}
