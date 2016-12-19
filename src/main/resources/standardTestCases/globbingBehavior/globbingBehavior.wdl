task writeFiles {
  String dir
  command <<<
    mkdir ${dir}
    echo "hullo" > ${dir}/hello.txt
    echo "buh-bye" > ${dir}/ciao.txt

    echo "what" > what.txt
    echo "how" > how.txt
    echo "why" > why.txt

    sleep 2
  >>>
  output {
    Array[File] strInterpolation=  glob("${dir}/*.txt")
    Array[File] callDirectory = glob("*.txt")
    Array[File] specificDirectory = glob(dir + "/*.txt")
  }
  runtime {
    docker:"ubuntu:latest"
  }
}

task unexpectedReturnCode {

  command {
    echo "easiest glob job ever!" > tiny.txt

    sleep 2
    exit 3
  }

  output {
    Array[File] easyGlob = glob("*.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    continueOnReturnCode: 3
  }
}

workflow globbingBehavior {
  String outputDir = "out"

  call writeFiles { input: dir = outputDir }

  call unexpectedReturnCode {}

  scatter (f in writeFiles.strInterpolation) {
    String strInterpolation_contents = read_string(f)
  }

  scatter (f in writeFiles.callDirectory) {
    String callDirectory_contents = read_string(f)
  }

  scatter (f in writeFiles.specificDirectory) {
    String specificDirectry_contents = read_string(f)
  }

  scatter (f in unexpectedReturnCode.easyGlob) {
    String easyGlob_contents = read_string(f)
  }

 output {
   strInterpolation_contents
   callDirectory_contents
   specificDirectory_contents
   easyGlob_contents
 }
}
