task validate_string {
    command {
        exit 0
    }
    output {
        String string_output = "hello world"
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

task validate_int {
    command {
        exit 0
    }
    output {
        Int int_output = 5
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

task validate_boolean {
    command {
        exit 0
    }
    output {
        Boolean boolean_output = true
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

task validate_float {
    command {
        exit 0
    }
    output {
        Float float_output = 5.5
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

workflow metadata_type_validation {
    call validate_string
    call validate_int
    call validate_boolean
    call validate_float
    output {
        validate_string.string_output
        validate_int.int_output
        validate_boolean.boolean_output
        validate_float.float_output
    }
}