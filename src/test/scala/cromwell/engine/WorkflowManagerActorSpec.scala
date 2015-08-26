package cromwell.engine

import java.util.UUID

import akka.testkit.TestActorRef
import cromwell.binding._
import cromwell.binding.command.CommandPart
import cromwell.binding.types.WdlStringType
import cromwell.binding.values.WdlString
import cromwell.engine.ExecutionStatus.{NotStarted, Running}
import cromwell.engine.backend.StdoutStderr
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.DataAccess.{WorkflowInfo, _}
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.{CallOutputs, WorkflowOutputs, _}
import cromwell.parser.BackendType
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.{HelloWorld, HelloWorldWithoutWorkflow, Incr}
import cromwell.{CromwellSpec, CromwellTestkitSpec, binding}

import scala.concurrent.duration.{Duration, _}
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

class WorkflowManagerActorSpec extends CromwellTestkitSpec("WorkflowManagerActorSpec") {

  "A WorkflowManagerActor" should {

    "run the Hello World workflow" in {
      withDataAccess { dataAccess =>
        implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance), self, "Test the WorkflowManagerActor")

        val workflowId = waitForHandledMessagePattern(pattern = "transitioning from Running to Succeeded") {
          messageAndWait[WorkflowId](SubmitWorkflow(HelloWorld.wdlSource(), HelloWorld.wdlJson, HelloWorld.rawInputs))
        }

        val status = messageAndWait[Option[WorkflowState]](WorkflowStatus(workflowId)).get
        status shouldEqual WorkflowSucceeded

        val workflowOutputs = messageAndWait[binding.WorkflowOutputs](WorkflowOutputs(workflowId))

        val actualWorkflowOutputs = workflowOutputs.map { case (k, WdlString(string)) => k -> string }
        actualWorkflowOutputs shouldEqual Map(HelloWorld.OutputKey -> HelloWorld.OutputValue)

        val callOutputs = messageAndWait[binding.CallOutputs](CallOutputs(workflowId, "hello.hello"))
        val actualCallOutputs = callOutputs.map { case (k, WdlString(string)) => k -> string }
        actualCallOutputs shouldEqual Map("salutation" -> HelloWorld.OutputValue)
      }
    }

    "Not try to restart any workflows when there are no workflows in restartable states" in {
      withDataAccess { dataAccess =>
        waitForPattern("Found no workflows to restart.") {
          TestActorRef(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance), self, "No workflows")
        }
      }
    }

    "Try to restart workflows when there are workflows in restartable states" in {
      val workflows = Map(
        WorkflowId(UUID.randomUUID()) -> WorkflowSubmitted,
        WorkflowId(UUID.randomUUID()) -> WorkflowRunning)
      val ids = workflows.keys.map(_.toString).toSeq.sorted
      val key = SymbolStoreKey("hello.hello", "addressee", None, input = true)
      val symbols = Map(key -> new SymbolStoreEntry(key, WdlStringType, Option(WdlString("world"))))

      withDataAccess { dataAccess =>
        import ExecutionContext.Implicits.global
        val setupFuture = Future.sequence(
          workflows map { case (workflowId, workflowState) =>
            val wdlSource = SampleWdl.HelloWorld.wdlSource()
            val wdlInputs = SampleWdl.HelloWorld.wdlJson
            val status = if (workflowState == WorkflowSubmitted) NotStarted else Running
            val workflowInfo = new WorkflowInfo(workflowId, wdlSource, wdlInputs)
            // FIXME? null AST
            val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
            val call = new Call(None, key.scope, task, Set.empty[FullyQualifiedName], Map.empty, None)
            for {
              _ <- dataAccess.createWorkflow(workflowInfo, symbols.values, Seq(call), new LocalBackend())
              _ <- dataAccess.updateWorkflowState(workflowId, workflowState)
              _ <- dataAccess.setStatus(workflowId, Seq(ExecutionDatabaseKey(call.fullyQualifiedName, None)), status)
            } yield ()
          }
        )
        Await.result(setupFuture, Duration.Inf)

        waitForPattern("Restarting workflow IDs: " + ids.mkString(", ")) {
          waitForPattern("Found 2 workflows to restart.") {
            // Workflows are always set back to Submitted on restart.
            waitForPattern("transitioning from Submitted to Running.", occurrences = 2) {
              // Both the previously in-flight call and the never-started call should get started.
              waitForPattern("starting calls: hello.hello", occurrences = 2) {
                waitForPattern("transitioning from Running to Succeeded", occurrences = 2) {
                  TestActorRef(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance), self, "2 restartable workflows")
                }
              }
            }
          }
        }
      }
    }

    val TestExecutionTimeout = 5000 milliseconds

    "Handle coercion failures gracefully" in {
      withDataAccess { dataAccess =>
        within(TestExecutionTimeout) {
          implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance), self, "Test WorkflowManagerActor coercion failures")
          waitForErrorWithException("Workflow failed submission") {
            Try {
              messageAndWait[WorkflowId](SubmitWorkflow(Incr.wdlSource(), Incr.wdlJson, Incr.rawInputs))
            } match {
              case Success(_) => fail("Expected submission to fail with uncoercable inputs")
              case Failure(e) =>
                e.getMessage shouldBe "The following errors occurred while processing your inputs:\n\nCould not coerce value for 'incr.incr.val' into: WdlIntegerType"
            }
          }
        }
      }
    }

    "error when running a workflowless WDL" in {
      withDataAccess { dataAccess =>
        implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance), self, "Test a workflowless submission")
        Try(messageAndWait[WorkflowId](SubmitWorkflow(HelloWorldWithoutWorkflow.wdlSource(),
          HelloWorldWithoutWorkflow.wdlJson, HelloWorldWithoutWorkflow.rawInputs))) match {
          case Success(_) => fail("Expected submission to fail due to no runnable workflows")
          case Failure(e) => e.getMessage shouldBe "Namespace does not have a local workflow to run"
        }
      }
    }

    "error when asked for outputs of a nonexistent workflow" in {
      withDataAccess { dataAccess =>
        within(TestExecutionTimeout) {
          implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance),
            self, "Test WorkflowManagerActor output lookup failure")
          val id = WorkflowId(UUID.randomUUID())
          Try {
            messageAndWait[binding.WorkflowOutputs](WorkflowOutputs(id))
          } match {
            case Success(_) => fail("Expected lookup to fail with unknown workflow")
            case Failure(e) => e.getMessage shouldBe s"Workflow '$id' not found"
          }
        }
      }
    }

    "error when asked for call logs of a nonexistent workflow" in {
      withDataAccess { dataAccess =>
        within(TestExecutionTimeout) {
          implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance),
            self, "Test WorkflowManagerActor call log lookup failure")
          val id = WorkflowId.randomId()
          Try {
            messageAndWait[StdoutStderr](CallStdoutStderr(id, "foo.bar"))
          } match {
            case Success(_) => fail("Expected lookup to fail with unknown workflow")
            case Failure(e) => e.getMessage shouldBe s"Workflow '$id' not found"
          }
        }
      }
    }

    "error when asked for logs of a nonexistent workflow" in {
      withDataAccess { dataAccess =>
        within(TestExecutionTimeout) {
          implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance),
            self, "Test WorkflowManagerActor log lookup failure")
          val id = WorkflowId.randomId()
          Try {
            messageAndWait[Map[LocallyQualifiedName, StdoutStderr]](WorkflowStdoutStderr(id))
          } match {
            case Success(_) => fail("Expected lookup to fail with unknown workflow")
            case Failure(e) => e.getMessage shouldBe s"Workflow '$id' not found"
          }
        }
      }
    }
  }
}
