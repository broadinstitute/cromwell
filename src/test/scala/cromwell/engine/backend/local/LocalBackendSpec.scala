package cromwell.engine.backend.local

import java.nio.file.Paths

import cromwell.CromwellTestkitSpec
import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend._
import cromwell.engine.workflow.BackendCallKey
import cromwell.util.SampleWdl
import org.specs2.mock.Mockito
import wdl4s.WdlSource
import wdl4s.values.WdlValue
class LocalBackendSpec extends CromwellTestkitSpec with Mockito {

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
    val backend = new LocalBackend(system)
    val backendCall = backend.bindCall(descriptor, BackendCallKey(call, None), Map.empty[String, WdlValue], abortRegistrationFunction = None)
    backendCall.execute map { _.result } map {
      case FailedExecution(e, _) => if (expectSuccess) fail("A call in a failOnStderr test which should have succeeded has failed ", e)
      case SuccessfulBackendCallExecution(_, _, _, _, _) => if (!expectSuccess) fail("A call in a failOnStderr test which should have failed has succeeded")
      case SuccessfulFinalCallExecution => fail("A FinalCall shouldn't have run in this failOnStderr")
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
