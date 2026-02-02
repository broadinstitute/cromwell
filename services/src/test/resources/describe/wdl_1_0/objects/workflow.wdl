version 1.0

workflow o {
  input {
    Object obj
    Object obj_optional = {"a": "is for apple"}
    MyStruct my_struct
    MyStruct my_struct_optional = {"myString": "zardoz", "myInt": 58}
    OtherStruct other
  }
  output {
    Object output_obj = {"b": "is for barnacle"}
    MyStruct my_struct_output = {"myString": "asdf", "myInt": 45}
  }
}

struct MyStruct {
  String myString
  Int myInt
}

struct OtherStruct {
  String myString
  Int myInt
}
