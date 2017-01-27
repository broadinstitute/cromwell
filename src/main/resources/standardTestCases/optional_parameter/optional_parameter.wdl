task select_first_with_optionals {
	String? supplied
	String? wfSupplied
	String? unsupplied
	String wfSuppliedWithDefault = select_first([wfSupplied, "HAPPY_BIRTHDAY_RUCHI"])
	String suppliedWithDefault = select_first([supplied, "HAPPY_BIRTHDAY_RUCHI"])
	String unsuppliedWithDefault = select_first([unsupplied, "HAPPY_BIRTHDAY_RUCHI"])

	command {
		echo "${wfSuppliedWithDefault}"
		echo "${suppliedWithDefault}"
		echo "${unsuppliedWithDefault}"
		sleep 2
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
		sleep 2
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
		sleep 2
	}

	output {
		Array[String] out = read_lines(stdout())
	}
}

workflow optional_parameter {
	call select_first_with_optionals { input: supplied = "SUPPLIED" }
	call interpolation_additions { input: supplied = "SUPPLIED", suppliedInt = 5}
	call default_interpolations { input: supplied = "SUPPLIED" }

	output {
	  Array[String] select_first_with_optionals_out = select_first_with_optionals.out
	  Int select_first_with_optionals_out_length = length(select_first_with_optionals_out)
	  Array[String] interpolation_additions_out = interpolation_additions.out
	  Int interpolation_additions_out_length = length(interpolation_additions_out)
	  Array[String] default_interpolations_out = default_interpolations.out
	  Int default_interpolations_out_length = length(default_interpolations_out)
	}
}
