version 1.0

workflow arrays {

  input {
    Array[Int] vanilla_array
    Array[Int]? optional_array
    Array[Int?] array_of_optionals
    Array[Int?]? optional_array_of_optionals
  }

  output {
    Array[Int] output_vanilla_array = []
    Array[Int]? output_optional_array = []
    Array[Int?] output_array_of_optionals = []
    Array[Int?]? output_optional_array_of_optionals = []
  }

}
