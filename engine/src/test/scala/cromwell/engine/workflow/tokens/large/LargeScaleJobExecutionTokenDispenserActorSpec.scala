package cromwell.engine.workflow.tokens.large

import akka.actor.{ActorRef, ActorSystem}
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

  var tokenDispenserUnderTest: ActorRef = _
  after {
    if (tokenDispenserUnderTest != null) system.stop(tokenDispenserUnderTest)
  }

  it should "limit two workflows to a max concurrency of 10 with no hog factor" ignore {
    val maxConcurrencyToTest = 10
    val hogFactor = 1 // ie, no hog factor
    val totalJobsPerWorkflow = maxConcurrencyToTest + 1
    val tokenType = JobExecutionTokenType(backendName, Some(maxConcurrencyToTest), hogFactor)

    tokenDispenserUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(maxConcurrencyToTest + 1, 100.millis)), "tokenDispenserUnderTest")

    val globalRunningJobsCounter = new RunningJobCounter()
    val bigWorkflow1 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), "multipleTokenUsingActor1")
    val bigWorkflow2 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), "multipleTokenUsingActor2")

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

    system.stop(bigWorkflow1)
    system.stop(bigWorkflow2)
  }

  it should "be able to restrain two workflows in the same hog group to a total of half of the total token pool" ignore {
    val totalTokensAvailable = 10
    val hogFactor = 2
    val maxConcurrencyExpected = 5
    val totalJobsPerWorkflow = maxConcurrencyExpected + 1
    val tokenType = JobExecutionTokenType(backendName, Some(totalTokensAvailable), hogFactor)

    tokenDispenserUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(maxConcurrencyExpected + 1, 100.millis)), "tokenDispenserUnderTest")

    val globalRunningJobsCounter = new RunningJobCounter()
    val bigWorkflow1 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), "multipleTokenUsingActor1")
    val bigWorkflow2 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), "multipleTokenUsingActor2")

    val parentProbe = new TestProbe(system, "parent")

    parentProbe.send(bigWorkflow1, Begin)
    parentProbe.send(bigWorkflow2, Begin)

    (0 until 2) foreach { _ =>
      parentProbe.expectMsgPF(100.seconds) {
        case TokenUsingActorCompletion(queueWaits, maximumConcurrency, errors) =>
          Assertions.assert(maximumConcurrency <= maxConcurrencyExpected, "(asserting maxActualConcurrency <= maxRequestedConcurrency)")
          queueWaits.size should be(totalJobsPerWorkflow)
          errors shouldBe List.empty
      }
    }

    globalRunningJobsCounter.getMax should be(maxConcurrencyExpected)
    system.stop(bigWorkflow1)
    system.stop(bigWorkflow2)
  }

  it should "be able to allocate two workflows in different hog groups exactly half of the total token pool each" ignore {
    val totalTokensAvailable = 10
    val hogFactor = 2
    val maxConcurrencyPerWorkflow = 5
    val maxConcurrencyOverall = 10
    val totalJobsPerWorkflow = maxConcurrencyPerWorkflow + 1
    val tokenType = JobExecutionTokenType(backendName, Some(totalTokensAvailable), hogFactor)

    tokenDispenserUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(maxConcurrencyOverall + 1, 100.millis)), "tokenDispenserUnderTest")

    val globalRunningJobsCounter = new RunningJobCounter()
    val bigWorkflow1 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), "multipleTokenUsingActor1")
    val bigWorkflow2 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupB", globalRunningJobsCounter), "multipleTokenUsingActor2")

    val parentProbe = new TestProbe(system, "parent")

    parentProbe.send(bigWorkflow1, Begin)
    parentProbe.send(bigWorkflow2, Begin)

    (0 until 2) foreach { _ =>
      parentProbe.expectMsgPF(100.seconds) {
        case TokenUsingActorCompletion(queueWaits, maximumConcurrency, errors) =>
          Assertions.assert(maximumConcurrency == maxConcurrencyPerWorkflow, "(asserting maxActualConcurrency <= maxRequestedConcurrency)")
          queueWaits.size should be(totalJobsPerWorkflow)
          errors shouldBe List.empty
      }
    }

    globalRunningJobsCounter.getMax should be(maxConcurrencyOverall)
    system.stop(bigWorkflow1)
    system.stop(bigWorkflow2)
  }

  it should "be able to allocate 100 workflows in 100 hog groups exactly 1/100 of the total token pool each" in {
    val totalTokensAvailable = 1000
    val hogFactor = 100
    val maxConcurrencyPerWorkflow = 10
    val maxConcurrencyOverall = 1000
    val totalJobsPerWorkflow = maxConcurrencyPerWorkflow * 2
    val tokenType = JobExecutionTokenType(backendName, Some(totalTokensAvailable), hogFactor)

    tokenDispenserUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(maxConcurrencyOverall + 1, 100.millis)), "tokenDispenserUnderTest")

    val globalRunningJobsCounter = new RunningJobCounter()

    val workflows = (0 until 100) map { i =>
      TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = s"hogGroup$i", globalRunningJobsCounter), s"multipleTokenUsingActor$i")
    }

    val parentProbe = new TestProbe(system, "parent")

    workflows foreach { parentProbe.send(_, Begin)}

    workflows.indices foreach { _ =>
      parentProbe.expectMsgPF(100.seconds) {
        case TokenUsingActorCompletion(queueWaits, maximumConcurrency, errors) =>
          Assertions.assert(maximumConcurrency == maxConcurrencyPerWorkflow, "(asserting maxActualConcurrency <= maxRequestedConcurrency)")
          queueWaits.size should be(totalJobsPerWorkflow)
          errors shouldBe List.empty
      }
    }

    globalRunningJobsCounter.getMax should be(maxConcurrencyOverall)
    workflows foreach { system.stop(_) }
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
