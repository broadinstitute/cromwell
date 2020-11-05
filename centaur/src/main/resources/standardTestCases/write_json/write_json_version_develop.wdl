version development

workflow write_json_version_develop {
    input {
        Array[Int] indices = [0, 1, 2]
    }
    scatter (i in indices) {
        call create_single_object { input: i=i }
    }

    output{
        Boolean boolean_out = read_boolean(write_json(false))
        Int int_out = read_int(write_json(1234))
        String string_out = read_string(write_json("Mr. Bean"))
        String array_out = read_string(write_json(indices))
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
        Object out = object {index: i, name: "mr_bean"}
	}
    runtime {
        docker: "ubuntu:latest"
    }
}
