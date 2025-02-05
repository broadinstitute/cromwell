task one {
  Int vertAxis
    command {
        echo ${vertAxis/2}
    }
    output {
        Float semiVertAxis = read_int(stdout())
    }
    runtime {
       docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:53a002b59dfcd43b4d15e97c1acbeae035ddd1b31a106659a312e9fe65f00afa"
    }
}

task two {
  Int horAxis
    command {
        echo ${horAxis/2}
    }
    output {
        Float semiHorAxis = read_int(stdout())
    }
    runtime {
       docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:53a002b59dfcd43b4d15e97c1acbeae035ddd1b31a106659a312e9fe65f00afa"
    }
}

task area{
   Float semiVertAxis
   Float semiHorAxis
   Float pi = 3.14159

   command {
   		echo ${semiHorAxis*semiVertAxis*pi}
   }
   output {
		Float ellipseArea = read_float(stdout())
   }
   runtime {
      docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
   }
}

workflow cacheBetweenWFNoCost {
      call one {
       input: vertAxis = 5
      }
      call two {
        input: horAxis = 6
      }
      call area {
        input: semiVertAxis = one.semiVertAxis, semiHorAxis = two.semiHorAxis
      }
    output {
        area.ellipseArea
    }
}