version 1.0

workflow read_write_json_roundtrip {
    input {
        Array[Int] indices = [0, 1, 2]
        Map[String, String] map_input = {"genre": "comedy", "movie": "mr. bean"}
        Pair[Int, String] pair_input = (1, "abc")
        Array[Map[String, String]] map_array = [{"artist": "maroon 5", "song": "memories"}]
    }

    scatter (i in indices) {
        Object object_array = object {index: i, name: "mr_bean"}
        Object object_with_array = object {index: i, samples: ["s1", "s2"]}
    }

    output{
        Boolean boolean_out = read_json(write_json(false))
        Int int_out = read_json(write_json(1234))
        Float float_out = read_json(write_json(123.456))
        String string_out = read_json(write_json("Mr. Bean"))
        Array[Int] array_out = read_json(write_json(indices))
        Object map_out = read_json(write_json(map_input))
        Object pair_out = read_json(write_json(pair_input))
        Object object_out = read_json(write_json(object_array[0]))
        Array[Object] array_object_out = read_json(write_json(object_array))
        Array[Object] array_object_with_array_out = read_json(write_json(object_with_array))
        Array[Map[String, String]] object_out_as_map = read_json(write_json(map_array))
    }
}
