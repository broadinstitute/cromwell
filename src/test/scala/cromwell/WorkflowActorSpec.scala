package cromwell

import akka.actor.ActorSystem
import akka.actor.SupervisorStrategy.Stop
import akka.testkit._
import com.typesafe.config.ConfigFactory
import cromwell.binding.WdlBinding
import cromwell.engine.WorkflowActor
import cromwell.engine.WorkflowActor._

import scala.concurrent.duration._
import scala.language.postfixOps

object WorkflowActorSpec {
  val config = """
                 |akka {
                 |  loggers = ["akka.testkit.TestEventListener"]
                 |  loglevel = "WARNING"
                 |}
               """.stripMargin

  val HelloWdl =
    """
      |task hello {
      |  command {
      |    echo "Hello ${addressee}"
      |  }
      |  output {
      |    String salutation = read_string("stdout")
      |  }
      |}
      |
      |workflow hello {
      |  call hello
      |}
    """.stripMargin
}

// Copying from http://doc.akka.io/docs/akka/snapshot/scala/testkit-example.html#testkit-example
class WorkflowActorSpec extends CromwellSpec(ActorSystem("WorkflowActorSpec", ConfigFactory.parseString(WorkflowActorSpec.config))) {

  import cromwell.binding.WdlImplicits._

  val helloBinding = WdlBinding.process(WorkflowActorSpec.HelloWdl)

  def buildWorkflowActor = TestActorRef[WorkflowActor]

  def constructWorkflowActor: TestActorRef[WorkflowActor] = {
    val workflowActor = buildWorkflowActor
    workflowActor ! Construct(helloBinding, Map(Addressee -> "world".toWdlValue))
    expectMsgPF() {
      case Constructed => ()
    }
    workflowActor
  }

  override def afterAll() {
    shutdown()
  }

  val Addressee = "hello.hello.addressee"
  val TestExecutionTimeout = 500 milliseconds

  "A WorkflowActor" should {

    "fail construction with missing inputs" in {
      within(TestExecutionTimeout) {
        buildWorkflowActor ! Construct(helloBinding, Map.empty)

        expectMsgPF() {
          case ConstructionFailed(unsatisfiedInputs) =>
            unsatisfiedInputs should have size 1
            unsatisfiedInputs should contain key Addressee
        }
      }
    }

    "fail construction with inputs of the wrong types" in {
      within(TestExecutionTimeout) {
        buildWorkflowActor ! Construct(helloBinding, Map(Addressee -> 3.toWdlValue))

        expectMsgPF() {
          case ConstructionFailed(unsatisfiedInputs) =>
            unsatisfiedInputs should have size 1
            unsatisfiedInputs should contain key Addressee
        }
      }
    }

    "fail construction when already constructed" in {
      within(TestExecutionTimeout) {
        constructWorkflowActor ! Construct(helloBinding, Map(Addressee -> "world".toWdlValue))
        expectMsgPF() {
          case InvalidOperation(_) => ()
        }
      }
    }

    "construct properly with proper inputs" in {
      within(TestExecutionTimeout) {
        constructWorkflowActor
      }
    }

    "fail to start when not first constructed" in {
      within(TestExecutionTimeout) {
        buildWorkflowActor ! Start
        expectMsgPF() {
          case InvalidOperation(_) => ()
        }
      }
    }

    "start" in {
      within(TestExecutionTimeout) {
        val workflowActor = constructWorkflowActor
        workflowActor ! Start
        expectMsgPF() {
          case Started => ()
        }
      }
    }

    "fail to stop when not first constructed" in {
      within(TestExecutionTimeout) {
        buildWorkflowActor ! Stop
        expectMsgPF() {
          case InvalidOperation(_) => ()
        }
      }
    }

    "stop" in {
      within(TestExecutionTimeout) {
        val probe = TestProbe()
        val workflowActor = constructWorkflowActor
        probe watch workflowActor
        workflowActor ! Stop
        expectMsgPF() {
          case Stopped => ()
        }
        probe expectTerminated workflowActor
      }
    }
  }
}
