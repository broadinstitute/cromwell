package cromwell.engine.backend.local

import java.util.UUID

import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding.command.{StringCommandPart, Command}
import cromwell.binding.values.WdlValue
import cromwell.binding.{Task, Call, WorkflowDescriptor}
import org.scalatest.{FlatSpec, Matchers}
import org.scalatest.mock.MockitoSugar
import org.mockito.Mockito._


import scala.util.{Try, Failure, Success}

class LocalBackendSpec extends FlatSpec with Matchers with MockitoSugar {

  val stdoutCommand = "echo beach"
  val stderrCommand = "echo beach >&2"

  // Mock setup:
  val mockWorkflowDescriptor: WorkflowDescriptor = mock[WorkflowDescriptor]
  when(mockWorkflowDescriptor.name) thenReturn "Mocky"
  when(mockWorkflowDescriptor.id) thenReturn UUID.randomUUID
  val mockScopedLookupFunction: ScopedLookupFunction = mock[ScopedLookupFunction]

  def testFailOnStderr(failOnStderr: Boolean, command: String, expectSuccess: Boolean): Unit = {
    val call = mock[Call]
    val task = mock[Task]
    when(call.docker) thenReturn None
    when(call.failOnStderr) thenReturn failOnStderr
    when(call.task) thenReturn task
    when(task.outputs) thenReturn Seq()
    new LocalBackend().executeCommand(command, mockWorkflowDescriptor, call, mockScopedLookupFunction) match {
      case Failure(e) => if (expectSuccess) fail("A call in a failOnStderr test which should have succeeded has failed ", e)
      case Success(_) => if (!expectSuccess) fail("A call in a failOnStderr test which should have failed has succeeded")
    }
  }

  "LocalBackend" should "allow stdout if failOnStderr is set" in {
    testFailOnStderr(failOnStderr = true, stdoutCommand, expectSuccess = true)
  }

  "LocalBackend" should "not allow stderr if failOnStderr is set" in {
      testFailOnStderr(failOnStderr = true, stderrCommand, expectSuccess = false)
  }

  "LocalBackend" should "allow stdout if failOnStderr is not set" in {
    testFailOnStderr(failOnStderr = false, stdoutCommand, expectSuccess = true)
  }

  "LocalBackend" should "allow stderr if failOnStderr is not set" in {
    testFailOnStderr(failOnStderr = false, stderrCommand, expectSuccess = true)
  }
}
