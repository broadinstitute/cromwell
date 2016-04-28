package cromwell.backend.impl.local

import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}

object TestWorkflows {

  case class TestWorkflow(workflowDescriptor: BackendWorkflowDescriptor,
                          config: BackendConfigurationDescriptor,
                          expectedResponse: BackendJobExecutionResponse)

  val HelloWorld =
    """
      |task hello {
      |  String addressee = "you "
      |  command {
      |    echo "Hello ${addressee}!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |
      |  RUNTIME
      |}
      |
      |workflow hello {
      |  call hello
      |}
    """.stripMargin

  val GoodbyeWorld =
    """
      |task goodbye {
      |  command {
      |    sh -c "exit 1"
      |  }
      |  output {
      |    String out = read_string(stdout())
      |  }
      |}
      |
      |workflow goodbye {
      |  call goodbye
      |}
    """.stripMargin

  val InputFiles =
    """
      |task localize {
      |  File inputFileFromJson
      |  File inputFileFromCallInputs
      |  command {
      |    cat ${inputFileFromJson}
      |    echo ""
      |    cat ${inputFileFromCallInputs}
      |  }
      |  output {
      |    Array[String] out = read_lines(stdout())
      |  }
      |
      |  RUNTIME
      |}
      |
      |workflow localize {
      |  File workflowFile
      |  call localize { input: inputFileFromCallInputs = workflowFile }
      |}
    """.stripMargin

  val InputExpressions =
    """
      |task expression {
      |  Int intNumber
      |  Float floatNumber
      |  Float sum = 15 + floatNumber
      |  command {
      |    echo ${intNumber} > ints
      |    echo ${floatNumber} > floats
      |    echo ${sum} >> floats
      |  }
      |  output {
      |    Array[Int] outInts = read_lines("ints")
      |    Array[Float] outFloats = read_lines("floats")
      |  }
      |}
      |
      |workflow expression {
      |  String str = "31380"
      |  Float f1 = 26
      |  call expression { input: intNumber = str,
      |                         floatNumber = f1 + 22}
      |}
    """.stripMargin

  val Sleep10 =
    """
      |task abort {
      |  command {
      |    sleep 10
      |    echo "something after sleep"
      |  }
      |}
      |
      |workflow abort {
      |  call abort
      |}
    """.stripMargin

  val Scatter =
    """
      |task scattering {
      |  Int intNumber
      |  command {
      |    echo ${intNumber}
      |  }
      |  output {
      |    Int out = read_string(stdout())
      |  }
      |}
      |
      |workflow scattering {
      |  Array[Int] numbers = [1, 2, 3]
      |  scatter (i in numbers) {
      |     call scattering { input: intNumber = i }
      |  }
      |}
    """.stripMargin

  val EngineFunctions =
    """
      |task scattering {
      |  Int intNumber
      |  command {
      |    echo ${intNumber}
      |  }
      |  output {
      |    Int out = read_string(stdout())
      |    String outStr= read_string(stdout())
      |  }
      |}
      |
      |workflow scattering {
      |  Array[Int] numbers = [1, 2, 3]
      |  scatter (i in numbers) {
      |     call scattering { input: intNumber = i }
      |  }
      |}
    """.stripMargin
}
