task one {
  Int radius = 62
    command {
        echo ${radius*radius}
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
   }
   output {
		Float area = read_float(stdout())
		Float piCopy = pi
   }
   runtime {
      docker: "ubuntu:latest"
   }

}

task three{
	Int rad
	   command {
		 	echo ${rad*rad}
		}
		output {
			Int rSquared = read_int(stdout())
			Int rCopy = rad
		}
		runtime {
			docker: "ubuntu:latest"
		}
}

workflow cacheWithinWF {
   call one {
   }
   call two {
   	input: r2 = one.rSquared
   }
   call three {
    input: rad = one.rCopy
   }
   call two as four {
   	input: r2 = three.rSquared
   }

   output {
        two.area
        four.area
   }

}