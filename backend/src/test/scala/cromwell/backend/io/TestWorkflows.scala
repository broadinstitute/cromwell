package cromwell.backend.io

import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}

object TestWorkflows {

  case class TestWorkflow(workflowDescriptor: BackendWorkflowDescriptor,
                          config: BackendConfigurationDescriptor,
                          expectedResponse: BackendJobExecutionResponse)

  val HelloWorld =
    s"""
      |task hello {
      |  String addressee = "you "
      |  command {
      |    echo "Hello $${addressee}!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |
      |  RUNTIME
      |}
      |
      |workflow wf_hello {
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
      |workflow wf_goodbye {
      |  call goodbye
      |}
    """.stripMargin

  val InputFiles =
    s"""
      |task localize {
      |  File inputFileFromJson
      |  File inputFileFromCallInputs
      |  command {
      |    cat $${inputFileFromJson}
      |    echo ""
      |    cat $${inputFileFromCallInputs}
      |  }
      |  output {
      |    Array[String] out = read_lines(stdout())
      |  }
      |
      |  RUNTIME
      |}
      |
      |workflow wf_localize {
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
      |workflow wf_abort {
      |  call abort
      |}
    """.stripMargin

  val Scatter =
    s"""
      |task scattering {
      |  Int intNumber
      |  command {
      |    echo $${intNumber}
      |  }
      |  output {
      |    Int out = read_string(stdout())
      |  }
      |}
      |
      |workflow wf_scattering {
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
      |workflow wf_localize {
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
      |workflow wf_localize {
      |  call localize
      |}
    """.stripMargin
  }
}
