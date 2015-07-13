package cromwell

import java.util.UUID

import akka.actor.{ActorRef, ActorSystem}
import akka.pattern.ask
import akka.testkit._
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec._
import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.DataAccess
import cromwell.engine.workflow.WorkflowActor
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.{WorkflowRunning, WorkflowSubmitted, WorkflowSucceeded}
import cromwell.util.SampleWdl
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.reflect.ClassTag

object CromwellTestkitSpec {
  val ConfigText =
    """
      |akka {
      |  loggers = ["akka.event.slf4j.Slf4jLogger", "akka.testkit.TestEventListener"]
      |  loglevel = "DEBUG"
      |  actor {
      |    debug {
      |       receive = on
      |    }
      |  }
      |}
    """.stripMargin

  implicit val timeout = Timeout(5 seconds)
}

abstract class CromwellTestkitSpec(name: String) extends TestKit(ActorSystem(name, ConfigFactory.parseString(ConfigText)))
with DefaultTimeout with ImplicitSender with WordSpecLike with Matchers with BeforeAndAfterAll with ScalaFutures {

  def startingCallsFilter[T](callNames: String*)(block: => T): T =
    waitForPattern(s"starting calls: ${callNames.mkString(", ")}$$")(block)

  def waitForHandledMessage[T](named: String)(block: => T): T = {
    waitForHandledMessagePattern(s"^received handled message $named")(block)
  }

  def waitForHandledMessagePattern[T](pattern: String)(block: => T): T = {
    EventFilter.debug(pattern = pattern, occurrences = 1).intercept {
      block
    }
  }

  val dataAccess = DataAccess()

  override protected def afterAll() = {
    super.afterAll()
    dataAccess.shutdown()
  }

  /**
   * Performs the following steps:
   *
   * <ol>
   * <li> Sends the specified message to the implicitly passed `ActorRef` via an `ask`.
   * <li> Collects the `Future[Any]` response.
   * <li> Downcasts the `Future[Any]` to a `Future[M]`.
   * <li> Issues a blocking `Await.result` on the `Future`, yielding an `M`.
   * </ol>
   *
   */
  def messageAndWait[M: ClassTag](message: AnyRef)(implicit actorRef: ActorRef): M = {
    val futureAny = actorRef ? message
    Await.result(futureAny.mapTo[M], 5 seconds)
  }

  /**
   * Wait for exactly one occurrence of the specified pattern in the specified block.  The block
   * is in its own parameter list for usage syntax reasons.
   */
  def waitForPattern[T](pattern: String, occurrences: Int = 1)(block: => T): T = {
    EventFilter.info(pattern = pattern, occurrences = occurrences).intercept {
      block
    }
  }

  private def buildFsmWorkflowActor(sampleWdl: SampleWdl, runtime: String) = {
    val namespace = WdlNamespace.load(sampleWdl.wdlSource(runtime))
    // This is a test and is okay with just throwing if coerceRawInputs returns a Failure.
    val coercedInputs = namespace.coerceRawInputs(sampleWdl.rawInputs).get
    val descriptor = WorkflowDescriptor(UUID.randomUUID(), namespace, sampleWdl.wdlSource(runtime), sampleWdl.wdlJson, coercedInputs)
    TestFSMRef(new WorkflowActor(descriptor, new LocalBackend, dataAccess))
  }

  def runWdlAndAssertOutputs(sampleWdl: SampleWdl, event: EventFilter, runtime: String = "", expectedOutputs: Map[FullyQualifiedName, WdlValue] = Map.empty): Unit = {
    val fsm = buildFsmWorkflowActor(sampleWdl, runtime)
    assert(fsm.stateName == WorkflowSubmitted)

    event.intercept {
      fsm ! Start
      within(5 seconds) {
        awaitCond(fsm.stateName == WorkflowRunning)
        awaitCond(fsm.stateName == WorkflowSucceeded)
        val outputs = fsm.ask(GetOutputs).mapTo[WorkflowOutputs].futureValue

        expectedOutputs foreach { case (outputFqn, expectedValue) =>
          val actualValue = outputs.getOrElse(outputFqn, throw new RuntimeException(s"Output $outputFqn not found"))
          actualValue shouldEqual expectedValue
        }
      }
    }
  }
}
