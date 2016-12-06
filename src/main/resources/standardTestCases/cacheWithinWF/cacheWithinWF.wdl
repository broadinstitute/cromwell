task one {
  Int radius
    command {
        echo ${radius*radius}
        sleep 2
    }
    output {
        Int rSquared = read_int(stdout())
		Int rCopy = radius
    }
    runtime {
       docker: "ubuntu:latest"
    }
}

task two{
   Int r2
   Float pi = 3.14159

   command {
   		echo ${r2*pi}
   		sleep 2
   }
   output {
		Float area = read_float(stdout())
		Float piCopy = pi
		Int rSquaredCopy = r2
   }
   runtime {
      docker: "ubuntu:latest"
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
