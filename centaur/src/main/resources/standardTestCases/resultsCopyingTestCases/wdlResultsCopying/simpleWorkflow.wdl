version 1.0

workflow simpleWorkflow {
    call simpleStdoutTask {}
    output {
        File outFile = simpleStdoutTask.outFile
        Array[File] outGlob = simpleStdoutTask.outGlob
        FooStruct myFoo = object {
            simple: 5,
            complex: ([42], { "t": simpleStdoutTask.objectOutputFile })
        }
        Array[Array[File]] nested_array = [[simpleStdoutTask.arraytest1, simpleStdoutTask.arraytest2], [simpleStdoutTask.arraytest3]]
    }
}

struct FooStruct {
  Int simple
  Pair[Array[Int], Map[String, File]] complex
}

task simpleStdoutTask {
  String outputFileName = "output.txt"
  String objectOutputFileName = "object_out.txt"

  command {
    echo 'Hello world' > ${outputFileName}
    echo 'Hello world' > ${objectOutputFileName}
    echo 'foo' > "foo.zardoz"
    echo 'bar' > "bar.zardoz"
    echo 'baz' > "baz.zardoz"

    echo "arraytest1" > "arraytest1.txt"
    echo "arraytest2" > "arraytest2.txt"
    echo "arraytest3" > "arraytest3.txt"
  }

  runtime {
    docker: "ubuntu"
    memory: "2G"
    cpu: 1
  }

  output {
    File outFile = outputFileName
    File objectOutputFile = objectOutputFileName
    Array[File] outGlob = glob("*.zardoz")
    File arraytest1 = "arraytest1.txt"
    File arraytest2 = "arraytest2.txt"
    File arraytest3 = "arraytest3.txt"
  }
}
