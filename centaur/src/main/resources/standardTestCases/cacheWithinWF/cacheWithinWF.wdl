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
       docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:53a002b59dfcd43b4d15e97c1acbeae035ddd1b31a106659a312e9fe65f00afa"
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
      docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:53a002b59dfcd43b4d15e97c1acbeae035ddd1b31a106659a312e9fe65f00afa"
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
