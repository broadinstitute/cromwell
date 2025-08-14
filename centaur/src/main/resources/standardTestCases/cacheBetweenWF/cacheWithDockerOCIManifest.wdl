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
       docker: "us-central1-docker.pkg.dev/broad-dsde-cromwell-dev/public-test-images/cromwell-oci-manifest-test:latest"
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
       docker: "us-central1-docker.pkg.dev/broad-dsde-cromwell-dev/public-test-images/cromwell-oci-manifest-test:latest"
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
