package cromwell

import java.util.UUID

import akka.actor.ActorSystem
import akka.testkit._
import com.typesafe.config.ConfigFactory
import cromwell.HelloWorldActorSpec._
import cromwell.binding.values.{WdlInteger, WdlString, WdlValue}
import cromwell.binding.{FullyQualifiedName, WdlBinding}
import cromwell.engine.WorkflowActor._
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.{UnsatisfiedInputsException, WorkflowActor}

import scala.concurrent.duration._
import scala.language.postfixOps

object HelloWorldActorSpec {
  val Config =
    """
      |akka {
      |  loggers = ["akka.testkit.TestEventListener"]
      |  loglevel = "DEBUG"
      |  actor.debug.receive = on
      |}
    """.stripMargin

  val HelloWdl =
    """
      |task hello {
      |  command {
      |    echo "Hello ${addressee}!"
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

  val Addressee = "hello.hello.addressee"
  val HelloInputs: Map[FullyQualifiedName, WdlValue] = Map(Addressee -> WdlString("world"))
}


// Copying from http://doc.akka.io/docs/akka/snapshot/scala/testkit-example.html#testkit-example
class HelloWorldActorSpec extends CromwellSpec(ActorSystem("HelloWorldActorSpec", ConfigFactory.parseString(Config))) {

  val helloBinding = WdlBinding.process(HelloWdl)

  def buildWorkflowActor(name: String = UUID.randomUUID().toString,
                         inputs: Map[FullyQualifiedName, WdlValue] = HelloInputs): TestActorRef[WorkflowActor] = {
    val binding = WdlBinding.process(HelloWdl)
    val props = WorkflowActor.buildWorkflowActorProps(UUID.randomUUID(), binding, inputs)
    TestActorRef(props, "Workflow-" + name)
  }

  override def afterAll() {
    shutdown()
  }

  val TestExecutionTimeout = 500 milliseconds

  "A WorkflowActor" should {

    "start" in {
      within(TestExecutionTimeout) {
        val workflowActor = buildWorkflowActor("started")
        startingCallsFilter("hello").intercept {
          workflowActor ! Start(new LocalBackend)
          expectMsgPF() {
            case Started => ()
          }
          expectMsgPF() {
            case Failed(t) =>
              fail(t)
            case Done(symbolStore) =>
              val maybeOutput = symbolStore.getOutputByFullyQualifiedName("hello.hello.salutation")

              val symbolStoreEntry = maybeOutput.getOrElse(throw new RuntimeException("No symbol store entry found!"))
              val wdlValue = symbolStoreEntry.wdlValue.getOrElse(throw new RuntimeException("No workflow output found!"))
              val actualOutput = wdlValue.asInstanceOf[WdlString].value.trim
              actualOutput shouldEqual "Hello world!"
          }
        }
      }
    }

    "fail to construct with missing inputs" in {
      intercept[UnsatisfiedInputsException] {
        buildWorkflowActor(inputs = Map.empty)
      }
    }

    "fail to construct with inputs of the wrong type" in {
      intercept[UnsatisfiedInputsException] {
        buildWorkflowActor(inputs = Map(Addressee -> WdlInteger(3)))
      }
    }
  }
}
