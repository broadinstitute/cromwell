version development

struct Random {
    Int index
    String name
}

workflow write_json_version_develop {
    input {
        Array[Int] indices = [0, 1, 2]
        Map[String, String] map_input = {"genre": "comedy", "movie": "mr. bean"}
        Pair[Int, String] pair_input = (1, "abc")
    }
    scatter (i in indices) {
        call create_single_object { input: i=i }
    }

    output{
        Boolean boolean_out = read_boolean(write_json(false))
        Int int_out = read_int(write_json(1234))
        Float float_out = read_float(write_json(123.456))
        String string_out = read_string(write_json("Mr. Bean"))
        String array_out = read_string(write_json(indices))
        String map_out = read_string(write_json(map_input))
        String pair_out = read_string(write_json(pair_input))
        String object_out = read_string(write_json(create_single_object.out[0]))
        String array_object_out = read_string(write_json(create_single_object.out))
    }
}

task create_single_object {
	input {
        Int i
	}
	command {
        echo "Creating single object"
	}
	output {
        Random out = object {index: i, name: "mr_bean"}
	}
    runtime {
        docker: "ubuntu:latest"
    }
}
