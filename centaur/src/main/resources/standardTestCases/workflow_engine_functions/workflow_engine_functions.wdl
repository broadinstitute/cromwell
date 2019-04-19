version 1.0


workflow w {
  String lines1 = "hello"
  String lines2 = "world"
  String map1 = "foo\tbar"
  String map2 = "baz\tqux"
  Array[String] array_lines = [lines1, lines2]
  Array[String] array_map = [map1, map2]
  File array_file = write_lines(array_lines)
  File map_file = write_lines(array_map)

  Array[String] out_array = read_lines(array_file)
  Map[String, String] out_map = read_map(map_file)

  output {
    String hello = out_array[0]
    String world = out_array[1]

    String bar = out_map["foo"]
    String qux = out_map["baz"]
  }
}
