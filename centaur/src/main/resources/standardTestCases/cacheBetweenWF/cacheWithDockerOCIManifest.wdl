task getAverage {
  Int base1 = 9
  Int base2 = 13
    command {
        echo ${(base1*base2)/2}
    }
    output {
        Float average = read_float(stdout())
    }
    runtime {
       docker: "gcr.io/broad-dsde-cromwell-dev/cromwell-oci-manifest-test@sha256:bac1423e2ae5b22a85905abff553cc9f0a7569ef0c1ab0732c821ce5e435a602"
    }
}

task heightProduct{
   Float baseAverage
   Int height = 7

   command {
   		echo ${baseAverage*height}
   }
   output {
		Float trapezoidalArea = read_float(stdout())
   }
   runtime {
       docker: "gcr.io/broad-dsde-cromwell-dev/cromwell-oci-manifest-test@sha256:bac1423e2ae5b22a85905abff553cc9f0a7569ef0c1ab0732c821ce5e435a602"
   }
}

workflow cacheWithDockerOCIManifest {
   call getAverage {
   }
   call heightProduct {
      input: baseAverage = getAverage.average
   }
   output {
        heightProduct.trapezoidalArea
   }
}
