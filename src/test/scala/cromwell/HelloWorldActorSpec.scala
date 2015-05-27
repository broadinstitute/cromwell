package cromwell

import java.util.UUID

import akka.actor.ActorSystem
import akka.testkit._
import com.typesafe.config.ConfigFactory
import cromwell.HelloWorldActorSpec._
import cromwell.binding.values.WdlString
import cromwell.binding.{UnsatisfiedInputsException, WdlBinding}
import cromwell.engine.WorkflowActor
import cromwell.engine.WorkflowActor._
import cromwell.engine.backend.local.LocalBackend
import cromwell.util.SampleWdl.HelloWorld
import cromwell.util.SampleWdl.HelloWorld.Addressee

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
}


// Copying from http://doc.akka.io/docs/akka/snapshot/scala/testkit-example.html#testkit-example
class HelloWorldActorSpec extends CromwellSpec(ActorSystem("HelloWorldActorSpec", ConfigFactory.parseString(Config))) {

  def buildWorkflowActor(name: String = UUID.randomUUID().toString,
                         rawInputs: binding.WorkflowRawInputs = HelloWorld.RawInputs): TestActorRef[WorkflowActor] = {
    val binding = WdlBinding.process(HelloWorld.WdlSource)
    val coercedInputs = binding.coerceRawInputs(rawInputs).get
    val props = WorkflowActor.buildWorkflowActorProps(UUID.randomUUID(), binding, coercedInputs)
    TestActorRef(props, self, "Workflow-" + name)
  }

  override def afterAll() {
    shutdown()
  }

  val TestExecutionTimeout = 5000 milliseconds

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
        buildWorkflowActor(rawInputs = Map.empty)
      }
    }

    "fail to construct with inputs of the wrong type" in {
      intercept[UnsatisfiedInputsException] {
        buildWorkflowActor(rawInputs = Map(Addressee -> 3))
      }
    }
  }
}
