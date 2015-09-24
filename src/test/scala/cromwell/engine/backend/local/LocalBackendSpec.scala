package cromwell.engine.backend.local

import java.nio.file.Paths

import cromwell.CromwellTestkitSpec
import cromwell.binding.values.WdlValue
import cromwell.binding.{WdlSource, WorkflowDescriptor}
import cromwell.engine.AbortRegistrationFunction
import cromwell.engine.backend.{AbortedExecution, FailedExecution, SuccessfulExecution}
import cromwell.engine.workflow.CallKey
import cromwell.util.SampleWdl

class LocalBackendSpec extends CromwellTestkitSpec("LocalBackendSpec") {

  object StdoutWdl extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """task out { command { echo beach } RUNTIME }
        |workflow wf { call out }
      """.stripMargin.replaceAll("RUNTIME", runtime)
    override val rawInputs =  Map.empty[String, String]
  }

  object StderrWdl extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """task err { command { echo beach >&2 } RUNTIME }
        |workflow wf { call err }
      """.stripMargin.replaceAll("RUNTIME", runtime)
    override val rawInputs =  Map.empty[String, String]
  }

  def stdoutDescriptor(runtime: String) = buildWorkflowDescriptor(StdoutWdl, runtime)
  def stderrDescriptor(runtime: String) = buildWorkflowDescriptor(StderrWdl, runtime)

  val prefix = "this_is_a_prefix"
  val cwd = Paths.get(".")

  def testFailOnStderr(descriptor: WorkflowDescriptor, expectSuccess: Boolean): Unit = {
    val call = descriptor.namespace.workflow.calls.head
    // NOTE: This test currently does not require storage. If that changes, use a withDataAccess() lender.
    val backend = new LocalBackend(dataAccess = null)
    val backendCall = backend.bindCall(descriptor, CallKey(call, None, None), Map.empty[String, WdlValue], AbortRegistrationFunction(_ => ()))
    backendCall.execute match {
      case FailedExecution(e, _) => if (expectSuccess) fail("A call in a failOnStderr test which should have succeeded has failed ", e)
      case SuccessfulExecution(_) => if (!expectSuccess) fail("A call in a failOnStderr test which should have failed has succeeded")
      case AbortedExecution => fail("Not expecting this at all")
    }
  }

  "LocalBackend" should {
    "allow stdout if failOnStderr is set" in {
      testFailOnStderr(stdoutDescriptor("runtime {failOnStderr: true}"), expectSuccess = true)
    }

    "not allow stderr if failOnStderr is set" in {
      testFailOnStderr(stderrDescriptor("runtime {failOnStderr: true}"), expectSuccess = false)
    }

    "allow stdout if failOnStderr is not set" in {
      testFailOnStderr(stdoutDescriptor("runtime {failOnStderr: false}"), expectSuccess = true)
    }

    "allow stderr if failOnStderr is not set" in {
      testFailOnStderr(stderrDescriptor("runtime {failOnStderr: false}"), expectSuccess = true)
    }
  }
}
