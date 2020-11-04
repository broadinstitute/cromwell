version 1.0

workflow write_json {
    input {
        Array[Int] indices = [0, 1, 2]
    }
    scatter (i in indices) {
        call create_single_object { input: i=i }
    }

    call write_array {input: array=indices}
    call write_object {input: obj=create_single_object.out[0]}
    call write_array_objects {input: arrayOfObj=create_single_object.out}

    output{
        String array_out = read_string(write_array.out)
        String object_out = read_string(write_object.out)
        String array_object_out = read_string(write_array_objects.out)
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

task write_array {
    input {
        Array[Int] array
    }
    command {
        cat ~{write_json(array)}
    }
    output {
        String out = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

task write_object {
    input {
        Object obj
    }
    command {
        cat ~{write_json(obj)}
    }
    output {
        String out = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

task write_array_objects {
    input {
        Array[Object] arrayOfObj
    }
    command {
        cat ~{write_json(arrayOfObj)}
    }
    output {
        String out = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
    }
}
