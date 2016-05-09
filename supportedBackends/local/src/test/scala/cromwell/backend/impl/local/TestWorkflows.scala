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

  val OutputProcess = {
    """
      |task localize {
      |  File inputFile
      |  command {
      |    echo "Hello" > a
      |    mkdir dir
      |    echo "world" > dir/b
      |  }
      |  output {
      |    File o1 = "a"
      |    Array[File] o2 = ["a", "dir/b"]
      |    File o3 = inputFile
      |  }
      |}
      |
      |workflow localize {
      |  call localize
      |}
    """.stripMargin
  }

  val MissingOutputProcess = {
    """
      |task localize {
      |  command {
      |  }
      |  output {
      |    File o1 = "c"
      |  }
      |}
      |
      |workflow localize {
      |  call localize
      |}
    """.stripMargin
  }
}
