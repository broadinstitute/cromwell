task oneAgain {
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

task twoAgain{
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

task threeAgain {
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

workflow cacheBetweenWF {

   call oneAgain {
   }
   call twoAgain {
   	input: r2 = oneAgain.rSquared
   }
   call threeAgain {
    input: rad = oneAgain.rCopy
   }
   call twoAgain as fourAgain {
   	input: r2 = threeAgain.rSquared
   }

   output {
        twoAgain.area
        fourAgain.area
   }

}