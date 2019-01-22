version 1.0

workflow o {
  input {
    Object obj
    Object obj_optional = {"a": "is for apple"}
    MyStruct ms
    MyStruct ms_optional = {"myString": "zardoz", "myInt": 58}
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
