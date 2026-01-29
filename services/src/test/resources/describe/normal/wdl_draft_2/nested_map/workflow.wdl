workflow wdl_draft_2_map {

  Map[Int, Map[Map[String, Int], String]] input_1

  output {
     Map[Int, Map[Map[String, Int], String]] output_1 = input_1
  }

}
