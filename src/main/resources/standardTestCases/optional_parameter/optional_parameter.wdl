task select_first_with_optionals {
	String? supplied
	String? wfSupplied
	String? unsupplied
	String suppliedWithDefault = select_first([supplied, "HAPPY_BIRTHDAY_RUCHI"])
	String wfSuppliedWithDefault = select_first([wfSupplied, "HAPPY_BIRTHDAY_RUCHI"])
	String unsuppliedWithDefault = select_first([unsupplied, "HAPPY_BIRTHDAY_RUCHI"])

	command {
		echo "${suppliedWithDefault}"
		echo "${wfSuppliedWithDefault}"
		echo "${unsuppliedWithDefault}"
	}
    runtime {
	  docker: "ubuntu:latest"
	}
	output {
		Array[String] out = read_lines(stdout())
	}
}

task optionals_with_defaults {
	String? supplied = "HAPPY_BIRTHDAY_RUCHI"
	String? call_and_wf_supplied = "HAPPY_BIRTHDAY_RUCHI"
	String? wfSupplied = "HAPPY_BIRTHDAY_RUCHI"
	String? unsupplied = "HAPPY_BIRTHDAY_RUCHI"

	command {
		echo "${supplied}"
		echo "${call_and_wf_supplied}"
		echo "${wfSupplied}"
		echo "${unsupplied}"
	}
    runtime {
	  docker: "ubuntu:latest"
	}
	output {
		Array[String] out = read_lines(stdout())
	}
}

task optional_array_interpolations {
	Array[Int]? supplied = [1, 2, 3]
	Array[Int]? defaulted = [1, 2, 3]
	Array[Int]? unsupplied

    String supplied_prefix = if defined(supplied) then "--foo " else ""
    String defaulted_prefix = if defined(defaulted) then "--foo " else ""
    String unsupplied_prefix = if defined(unsupplied) then "--foo " else ""

	command {
		echo "${supplied_prefix}${sep="," supplied}"
		echo "${defaulted_prefix}${sep="," defaulted}"
		echo "${unsupplied_prefix}${sep="," unsupplied}"
	}
    runtime {
	  docker: "ubuntu:latest"
	}
	output {
		Array[String] out = read_lines(stdout())
	}
}

task interpolation_additions {
	String? supplied
	Int? suppliedInt
	String? wfSupplied
	String? unsupplied
	String? unsuppliedInt

	command {
		echo ${"supplied: " + supplied}
		echo ${"suppliedInt: " + suppliedInt}
		echo ${"wfSupplied: " + wfSupplied}
		echo ${"unsupplied: " + unsupplied}
		echo ${"unsuppliedInt: " + unsuppliedInt}
	}
	runtime {
        docker: "ubuntu:latest"
    }
	output {
		Array[String] out = read_lines(stdout())
	}
}

task default_interpolations {
	String? supplied
	String? wfSupplied
	String? unsupplied

	command {
		echo ${default="HAPPY_BIRTHDAY_RUCHI" supplied}
		echo ${default="HAPPY_BIRTHDAY_RUCHI" wfSupplied}
		echo ${default="HAPPY_BIRTHDAY_RUCHI" unsupplied}
	}
	runtime {
        docker: "ubuntu:latest"
    }
	output {
		Array[String] out = read_lines(stdout())
	}
}

workflow optional_parameter {
	call select_first_with_optionals { input: supplied = "SUPPLIED" }
	call optionals_with_defaults { input: supplied = "SUPPLIED", call_and_wf_supplied = "SUPPLIED" }
	call optional_array_interpolations { input: supplied = [5, 5, 5] }
	call interpolation_additions { input: supplied = "SUPPLIED", suppliedInt = 5}
	call default_interpolations { input: supplied = "SUPPLIED" }

	output {
	  Array[String] select_first_with_optionals_out = select_first_with_optionals.out
	  Int select_first_with_optionals_out_length = length(select_first_with_optionals_out)

	  Array[String] optionals_with_defaults_out = optionals_with_defaults.out
	  Int optionals_with_defaults_out_length = length(optionals_with_defaults_out)

	  Array[String] optional_array_interpolations_out = optional_array_interpolations.out
      Int optional_array_interpolations_out_length = length(optional_array_interpolations_out)

	  Array[String] interpolation_additions_out = interpolation_additions.out
	  Int interpolation_additions_out_length = length(interpolation_additions_out)

	  Array[String] default_interpolations_out = default_interpolations.out
	  Int default_interpolations_out_length = length(default_interpolations_out)
	}
}
