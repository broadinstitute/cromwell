package cromwell.engine.backend.local

import java.io.File
import java.nio.file.{Files, Paths}

import akka.testkit._
import cromwell.CromwellTestkitSpec
import cromwell.binding.values.{WdlFile, WdlValue}
import cromwell.binding.{WdlSource, WorkflowDescriptor}
import cromwell.engine.AbortRegistrationFunction
import cromwell.engine.backend.Backend.StdoutStderrException
import cromwell.util.SampleWdl

import scala.util.{Failure, Success}

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
    val backend = new LocalBackend()
    val backendCall = backend.bindCall(descriptor, call, Map.empty[String, WdlValue], AbortRegistrationFunction(_ => ()))
    backendCall.execute match {
      case Failure(e) => if (expectSuccess) fail("A call in a failOnStderr test which should have succeeded has failed ", e)
      case Success(_) => if (!expectSuccess) fail("A call in a failOnStderr test which should have failed has succeeded")
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
