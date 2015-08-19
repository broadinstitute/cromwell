package cromwell.engine.backend.local

import java.io.File
import java.nio.file.{Files, Paths}
import java.util.UUID

import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding.values.WdlFile
import cromwell.binding.{Call, Task, WorkflowDescriptor}
import cromwell.engine.backend.Backend.StdoutStderrException
import cromwell.engine.{AbortFunctionRegistration, WorkflowId}
import org.mockito.Mockito._
import org.scalatest.mock.MockitoSugar
import org.scalatest.{FlatSpec, Matchers}

import scala.util.{Failure, Success}

class LocalBackendSpec extends FlatSpec with Matchers with MockitoSugar {

  val prefix = "this_is_a_prefix"
  val cwd = Paths.get(".")
  val stdoutCommand = "echo beach"
  val stderrCommand = "echo beach >&2"

  // Mock setup:
  val mockWorkflowDescriptor: WorkflowDescriptor = mock[WorkflowDescriptor]
  when(mockWorkflowDescriptor.name) thenReturn "Mocky"
  when(mockWorkflowDescriptor.id) thenReturn WorkflowId(UUID.randomUUID)
  val mockScopedLookupFunction: ScopedLookupFunction = mock[ScopedLookupFunction]

  def testFailOnStderr(failOnStderr: Boolean, command: String, expectSuccess: Boolean): Unit = {
    val call = mock[Call]
    val task = mock[Task]
    when(call.docker) thenReturn None
    when(call.failOnStderr) thenReturn failOnStderr
    when(call.task) thenReturn task
    when(task.outputs) thenReturn Seq()
    new LocalBackend().executeCommand(command, mockWorkflowDescriptor, call, call.inputMappings, mockScopedLookupFunction, AbortFunctionRegistration(_ => ())) match {
      case Failure(e) => if (expectSuccess) fail("A call in a failOnStderr test which should have succeeded has failed ", e)
      case Success(_) => if (!expectSuccess) fail("A call in a failOnStderr test which should have failed has succeeded")
    }
  }

  "LocalBackend" should "allow stdout if failOnStderr is set" in {
    testFailOnStderr(failOnStderr = true, stdoutCommand, expectSuccess = true)
  }

  it should "not allow stderr if failOnStderr is set" in {
      testFailOnStderr(failOnStderr = true, stderrCommand, expectSuccess = false)
  }

  it should "allow stdout if failOnStderr is not set" in {
    testFailOnStderr(failOnStderr = false, stdoutCommand, expectSuccess = true)
  }

  it should "allow stderr if failOnStderr is not set" in {
    testFailOnStderr(failOnStderr = false, stderrCommand, expectSuccess = true)
  }

  "LocalBackend temp file finder" should "be able to find a temp file" in {
    val fp = File.createTempFile(prefix, ".tmp", cwd.toFile).toPath
    val found = LocalBackend.findTempFile(cwd, "this_is_a_prefix")
    found shouldEqual WdlFile(fp.toAbsolutePath.toString)
    Files.delete(fp)
  }

  it should "return an error if more than one file match" in {
    val fp = File.createTempFile(prefix, ".tmp", cwd.toFile).toPath
    val fp2 = File.createTempFile(prefix, ".tmp", cwd.toFile).toPath
    try {
      LocalBackend.findTempFile(cwd, prefix)
      fail("Expected findTempFile to fail if more than one file matched")
    } catch {
      case e:StdoutStderrException => // expected
    }
    Files.delete(fp)
    Files.delete(fp2)
  }

  "LocalBackend temp file finder" should "return an error if no files match" in {
    try {
      LocalBackend.findTempFile(cwd, prefix)
      fail("Expected findTempFile to fail if no files matched")
    } catch {
      case e:StdoutStderrException => // expected
    }
  }
}
