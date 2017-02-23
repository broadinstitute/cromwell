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
       docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
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
      docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
   }
}

workflow cacheBetweenWF {
   call getAverage {
   }
   call heightProduct {
      input: baseAverage = getAverage.average
   }
   output {
        heightProduct.trapezoidalArea
   }
}
