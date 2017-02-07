task makeString {

	command {
		echo hello && sleep 2
	}
	runtime {
	  docker: "ubuntu:latest"
	}
	output {
	    String out = read_string(stdout())
	}
}


workflow selection_functions {
  String? created
  String? notCreated

  if (true) {
    call makeString as makeStringTrue
  }

  if (false) {
    call makeString as makeStringFalse
  }

  output {
    Array[String] allFromTruths = select_all([makeStringTrue.out, makeStringTrue.out])
    Array[String] allFromFalses = select_all([makeStringFalse.out, makeStringFalse.out])
    Int allFromFalsesLength = length(allFromFalses)
    Array[String] allFromMix = select_all([makeStringTrue.out, makeStringFalse.out])
    String firstFromTruths = select_first([makeStringTrue.out, makeStringTrue.out])
    String firstFromFalses = select_first([makeStringFalse.out, makeStringFalse.out, "default"])
    String firstFromMix = select_first([makeStringFalse.out, makeStringTrue.out])
    String? mst = makeStringTrue.out
    Array[String?] msts = [mst, mst]
    Array[String] allFromMsts = select_all(msts)
    String firstFromMsts = select_first(msts)
  }
}
