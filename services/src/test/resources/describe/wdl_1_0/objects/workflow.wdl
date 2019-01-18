version 1.0

workflow o {
  input {
    Object obj = {"a": "is for apple"}
    MyStruct ms
  }
  output {
    Object output_obj = {"b": "is for barnacle"}
    MyStruct ms_output = {"myString": "asdf", "myInt": 45}
  }
}

struct MyStruct {
  String myString
  Int myInt
}
