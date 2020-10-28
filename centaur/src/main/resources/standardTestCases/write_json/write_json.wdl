version 1.0

workflow write_json {
    input {
        Array[Int] indicies = [0, 1, 2]
    }
    scatter (i in indicies) {
        call create_single_object { input: i=i }
    }

    call write_array {input: array=indicies}
    call write_object {input: obj=create_single_object.out[0]}
    call write_array_objects {input: array=create_single_object.out}

    output{
        String array_out = write_array.out
        String object_out = write_object.out
        String array_object_out = write_array_objects.out
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
		Object out = object {name: "mr_bean", num: i}
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
        Array[Object] array
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
