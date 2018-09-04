package cromwell.engine.workflow.tokens.large

import akka.actor.ActorSystem
import akka.testkit.{ImplicitSender, TestActorRef, TestKit, TestProbe}
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.DynamicRateLimiter.Rate
import cromwell.engine.workflow.tokens.large.PatientTokenNeedingActor.Begin
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor
import cromwell.engine.workflow.tokens.large.LargeScaleJobExecutionTokenDispenserActorSpec.RunningJobCounter
import cromwell.engine.workflow.tokens.large.MultipleTokenUsingActor.TokenUsingActorCompletion
import org.scalatest._
import org.scalatest.concurrent.Eventually

import scala.concurrent.duration._

class LargeScaleJobExecutionTokenDispenserActorSpec extends TestKit(ActorSystem("LSJETDASpec")) with ImplicitSender with FlatSpecLike with Matchers with BeforeAndAfter with BeforeAndAfterAll with Eventually {

  val backendName = "PAPI"

  behavior of "JobExecutionTokenDispenserActor with concurrent demands"


  it should "correctly limit to a max concurrency of 10 with no hog factor" in {
    val maxConcurrencyToTest = 10
    val hogFactor = 1 // ie, no hog factor
    val totalJobsPerWorkflow = maxConcurrencyToTest + 1
    val tokenType = JobExecutionTokenType(backendName, Some(maxConcurrencyToTest), hogFactor)

    val actorRefUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(maxConcurrencyToTest + 1, 100.millis)), "tokenDispenserUnderTest")

    val globalRunningJobsCounter = new RunningJobCounter()
    val bigWorkflow1 = TestActorRef(new MultipleTokenUsingActor(actorRefUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), "multipleTokenUsingActor1")
    val bigWorkflow2 = TestActorRef(new MultipleTokenUsingActor(actorRefUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), "multipleTokenUsingActor2")

    val parentProbe = new TestProbe(system, "parent")

    parentProbe.send(bigWorkflow1, Begin)
    parentProbe.send(bigWorkflow2, Begin)

    (0 until 2) foreach { _ =>
      parentProbe.expectMsgPF(100.seconds) {
        case TokenUsingActorCompletion(queueWaits, maximumConcurrency, errors) =>
          Assertions.assert(maximumConcurrency <= maxConcurrencyToTest, "(asserting maxActualConcurrency <= maxRequestedConcurrency)")
          queueWaits.size should be(totalJobsPerWorkflow)
          errors shouldBe List.empty
      }
    }

    globalRunningJobsCounter.getMax should be(maxConcurrencyToTest)
  }

}

object LargeScaleJobExecutionTokenDispenserActorSpec {

  /**
    * Like an AtomicInteger, but it remembers the highest value it's ever seen
    */
  final class RunningJobCounter() {
    private var value: Int = 0
    private var max: Int = 0

    def increment(): Unit = {
      this.synchronized {
        value = value + 1
        max = Math.max(value, max)
      }
      ()
    }

    def decrement(): Unit = {
      this.synchronized {
        value = value - 1
      }
      ()
    }

    def getMax: Int = max
  }
}
