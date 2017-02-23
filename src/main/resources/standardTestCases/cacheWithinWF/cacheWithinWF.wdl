task one {
  Int radius
    command {
        echo ${radius*radius}
    }
    output {
        Int rSquared = read_int(stdout())
		Int rCopy = radius
    }
    runtime {
       docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    }
}

task two{
   Int r2
   Float pi = 3.14159

   command {
   		echo ${r2*pi}
   }
   output {
		Float area = read_float(stdout())
		Float piCopy = pi
		Int rSquaredCopy = r2
   }
   runtime {
      docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
   }
}

workflow cacheWithinWF {
   call one {
    input: radius = 62
   }
   call two {
   	input: r2 = one.rSquared
   }
   call two as twoAgain {
   	input: r2 = two.rSquaredCopy
   }
   output {
      two.area
      twoAgain.area
   }
}
