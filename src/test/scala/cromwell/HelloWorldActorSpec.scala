package cromwell

import java.util.UUID

import akka.actor.ActorSystem
import akka.testkit._
import akka.pattern.ask
import com.typesafe.config.ConfigFactory
import cromwell.HelloWorldActorSpec._
import cromwell.binding.values.WdlString
import cromwell.binding.{WorkflowOutputs, UnsatisfiedInputsException, WdlNamespace}
import cromwell.engine.{WorkflowSucceeded, WorkflowRunning, WorkflowSubmitted, WorkflowActor}
import cromwell.engine.WorkflowActor._
import cromwell.engine.backend.local.LocalBackend
import cromwell.util.SampleWdl.HelloWorld
import cromwell.util.SampleWdl.HelloWorld.Addressee

import scala.concurrent.Await
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
class HelloWorldActorSpec extends CromwellTestkitSpec(ActorSystem("HelloWorldActorSpec", ConfigFactory.parseString(Config))) {
  private def buildWorkflowActor(name: String = UUID.randomUUID().toString,
                         rawInputs: binding.WorkflowRawInputs = HelloWorld.RawInputs): TestActorRef[WorkflowActor] = {
    val namespace = WdlNamespace.load(HelloWorld.WdlSource)
    val coercedInputs = namespace.coerceRawInputs(rawInputs).get
    val props = WorkflowActor.props(UUID.randomUUID(), namespace, coercedInputs, new LocalBackend)
    TestActorRef(props, self, "Workflow-" + name)
  }

  override def afterAll() {
    shutdown()
  }

  val TestExecutionTimeout = 5000 milliseconds

  "A WorkflowActor" should {
    "start" in {
      val namespace = WdlNamespace.load(HelloWorld.WdlSource)
      val coercedInputs = namespace.coerceRawInputs(HelloWorld.RawInputs).get
      /*
        The TestFSMRef is kind of quirky, defining it here instead of the buildWorkflowActor function. It could
        be generalized a bit but it is probably not worth the hassle for a test class
       */
      val fsm = TestFSMRef(new WorkflowActor(UUID.randomUUID(), namespace, coercedInputs, new LocalBackend))
      assert(fsm.stateName == WorkflowSubmitted)
      startingCallsFilter("hello").intercept {
        fsm ! Start
        within(TestExecutionTimeout) {
          awaitCond(fsm.stateName == WorkflowRunning)
          awaitCond(fsm.stateName == WorkflowSucceeded)
          val outputName = "hello.hello.salutation"
          val outputs = Await.result(fsm.ask(GetOutputs).mapTo[WorkflowOutputs], 5 seconds)
          val salutation = outputs.getOrElse(outputName, throw new RuntimeException(s"Output '$outputName' not found."))
          val actualOutput = salutation.asInstanceOf[WdlString].value.trim
          actualOutput shouldEqual "Hello world!"
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
