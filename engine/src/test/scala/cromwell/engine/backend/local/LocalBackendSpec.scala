package cromwell.engine.backend.local

import java.nio.file.Paths

import cromwell.CromwellTestkitSpec
import cromwell.CromwellTestkitSpec.TestWorkflowManagerSystem
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.backend._
import cromwell.engine.workflow.BackendCallKey
import cromwell.util.SampleWdl
import org.scalatest.FlatSpec
import org.specs2.mock.Mockito
import wdl4s.WdlSource
import wdl4s.values.WdlValue

class LocalBackendSpec extends CromwellTestkitSpec with Mockito with WorkflowDescriptorBuilder {

  val testWorkflowManagerSystem = new TestWorkflowManagerSystem
  override implicit val actorSystem = testWorkflowManagerSystem.actorSystem

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

  def buildWorkflowWithRuntime(wdl: SampleWdl, runtime: String): OldStyleWorkflowDescriptor = {
    val sources = WorkflowSourceFiles(wdl.wdlSource(runtime), wdl.wdlJson, "{}")
    materializeWorkflowDescriptorFromSources(workflowSources = sources)
  }

  def stdoutDescriptor(runtime: String) = buildWorkflowWithRuntime(StdoutWdl, runtime)
  def stderrDescriptor(runtime: String) = buildWorkflowWithRuntime(StderrWdl, runtime)

  val prefix = "this_is_a_prefix"
  val cwd = Paths.get(".")

  def testFailOnStderr(descriptor: OldStyleWorkflowDescriptor, expectSuccess: Boolean): Unit = {
    val call = descriptor.namespace.workflow.calls.head
    val jobDescriptor = OldStyleBackendCallJobDescriptor(descriptor, BackendCallKey(call, None, 1), Map.empty[String, WdlValue])
    jobDescriptor.execute map { _.result } map {
      case OldStyleNonRetryableFailedExecution(e, _, _) => if (expectSuccess) fail("A call in a failOnStderr test which should have succeeded has failed ", e)
      case OldStyleRetryableFailedExecution(e, _, _) => if (expectSuccess) fail("A call in a failOnStderr test which should have succeeded has failed ", e)
      case OldStyleSuccessfulBackendCallExecution(_, _, _, _, _) => if (!expectSuccess) fail("A call in a failOnStderr test which should have failed has succeeded")
      case OldStyleSuccessfulFinalCallExecution => fail("A FinalCall shouldn't have run in this failOnStderr")
      case OldStyleAbortedExecution => fail("Not expecting this at all")
    }
  }

  "LocalBackend" should {
    "allow stdout if failOnStderr is set" ignore {
      testFailOnStderr(stdoutDescriptor("runtime {failOnStderr: true}"), expectSuccess = true)
    }

    "not allow stderr if failOnStderr is set" ignore {
      testFailOnStderr(stderrDescriptor("runtime {failOnStderr: true}"), expectSuccess = false)
    }

    "allow stdout if failOnStderr is not set" ignore {
      testFailOnStderr(stdoutDescriptor("runtime {failOnStderr: false}"), expectSuccess = true)
    }

    "allow stderr if failOnStderr is not set" ignore {
      testFailOnStderr(stderrDescriptor("runtime {failOnStderr: false}"), expectSuccess = true)
    }
  }
}
