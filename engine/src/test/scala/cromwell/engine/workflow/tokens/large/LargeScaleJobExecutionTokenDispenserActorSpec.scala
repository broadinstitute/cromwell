package cromwell.engine.workflow.tokens.large

import java.util.concurrent.atomic.AtomicInteger

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

  val multipleTokenUsingActorIndex: AtomicInteger = new AtomicInteger(0)
  def multipleTokenUsingActorName() = s"multipleTokenUsingActor${multipleTokenUsingActorIndex.getAndIncrement()}"

  val backendName = "PAPI"

  behavior of "JobExecutionTokenDispenserActor with concurrent demands"

  it should "limit two workflows to a max concurrency of 10 with no hog factor" in {
    val maxConcurrencyToTest = 10
    val hogFactor = 1 // ie, no hog factor
    val totalJobsPerWorkflow = maxConcurrencyToTest + 1
    val tokenType = JobExecutionTokenType(backendName, Some(maxConcurrencyToTest), hogFactor)

    val tokenDispenserUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(maxConcurrencyToTest + 1, 100.millis), None), "tokenDispenserUnderTest1")

    val globalRunningJobsCounter = new RunningJobCounter()
    val bigWorkflow1 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), multipleTokenUsingActorName())
    val bigWorkflow2 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), multipleTokenUsingActorName())

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
    system.stop(tokenDispenserUnderTest)
  }

  it should "be able to restrain two workflows in the same hog group to a total of half of the total token pool" in {
    val totalTokensAvailable = 10
    val hogFactor = 2
    val maxConcurrencyExpected = 5
    val totalJobsPerWorkflow = maxConcurrencyExpected + 1
    val tokenType = JobExecutionTokenType(backendName, Some(totalTokensAvailable), hogFactor)

    val tokenDispenserUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(maxConcurrencyExpected + 1, 100.millis), None), "tokenDispenserUnderTest2")

    val globalRunningJobsCounter = new RunningJobCounter()
    val bigWorkflow1 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), multipleTokenUsingActorName())
    val bigWorkflow2 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), multipleTokenUsingActorName())

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
    system.stop(tokenDispenserUnderTest)
  }

  it should "be able to allocate two workflows in two hog groups exactly half of the total token pool each" in {
    val totalTokensAvailable = 10
    val hogFactor = 2
    val maxConcurrencyPerWorkflow = 5
    val maxConcurrencyOverall = 10
    val totalJobsPerWorkflow = maxConcurrencyPerWorkflow + 1
    val tokenType = JobExecutionTokenType(backendName, Some(totalTokensAvailable), hogFactor)

    val tokenDispenserUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(maxConcurrencyOverall + 1, 100.millis), None), "tokenDispenserUnderTest3")

    val globalRunningJobsCounter = new RunningJobCounter()
    val bigWorkflow1 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupA", globalRunningJobsCounter), multipleTokenUsingActorName())
    val bigWorkflow2 = TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = "hogGroupB", globalRunningJobsCounter), multipleTokenUsingActorName())

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
    system.stop(tokenDispenserUnderTest)
  }

  it should "be able to allocate 100 workflows in 100 hog groups exactly 1/100 of the total token pool each" in {
    val totalTokensAvailable = 1000
    val hogFactor = 100
    val maxConcurrencyPerWorkflow = 10
    val maxConcurrencyOverall = 1000
    val totalJobsPerWorkflow = maxConcurrencyPerWorkflow * 2
    val tokenType = JobExecutionTokenType(backendName, Some(totalTokensAvailable), hogFactor)

    val tokenDispenserUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(maxConcurrencyOverall + 1, 100.millis), None), "tokenDispenserUnderTest4")

    val globalRunningJobsCounter = new RunningJobCounter()

    val workflows = (0 until 100) map { i =>
      TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = s"hogGroup$i", globalRunningJobsCounter), multipleTokenUsingActorName())
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
    system.stop(tokenDispenserUnderTest)
  }

  it should "be able to allocate 100 workflows in 25 hog groups exactly 1/50 of the total token pool per hog group" in {
    val totalTokensAvailable = 1000
    val totalWorkflows = 100
    val totalHogGroups = 25
    val hogFactor = 50
    val maxConcurrencyPerHogGroup = 20
    val maxConcurrencyOverall = 1000
    val totalJobsPerWorkflow = maxConcurrencyPerHogGroup * 2
    val tokenType = JobExecutionTokenType(backendName, Some(totalTokensAvailable), hogFactor)

    val tokenDispenserUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(maxConcurrencyOverall + 1, 100.millis), None), "tokenDispenserUnderTest5")

    val hogGroupConcurrencyCounters = (0 until totalHogGroups).toVector map { _ => new RunningJobCounter() }

    val workflows = (0 until totalWorkflows) map { i =>
      val hogGroupNumber = i % totalHogGroups
      TestActorRef(new MultipleTokenUsingActor(tokenDispenserUnderTest, tokenType, totalJobsPerWorkflow, hogGroup = s"hogGroup$hogGroupNumber", hogGroupConcurrencyCounters(hogGroupNumber)), multipleTokenUsingActorName())
    }

    val parentProbe = new TestProbe(system, "parent")

    workflows foreach { parentProbe.send(_, Begin)}

    workflows.indices foreach { _ =>
      parentProbe.expectMsgPF(100.seconds) {
        case TokenUsingActorCompletion(queueWaits, maximumConcurrency, errors) =>
          Assertions.assert(maximumConcurrency <= maxConcurrencyPerHogGroup, "(asserting maxActualConcurrency per workflow <= maxRequestedConcurrency per hog group)")
          queueWaits.size should be(totalJobsPerWorkflow)
          errors shouldBe List.empty
      }
    }

    // Check that each hog group was able to start jobs all the way up to its limit (but not beyond):
    // NOTE:
    // The test is set up to let the workflows of every hog group reach exactly the concurrency limit in parallel, but because we're in a non-deterministic test,
    // occasionally it gets off-by-one (eg one hog group maxing out at 39 instead of 40 concurrent jobs).
    // So:
    //   - check that at least 95% of the hog groups did exactly hit the limit (as we would expect):
    //   - check that no hog groups ever go over the limit,
    //   - check that all hog groups got to at least 95% of what we expected, even if they didn't quite reach the exact limit
    var exactlyAtLimit = 0
    (0 until totalHogGroups).toVector foreach { hogGroupNumber =>
      val c: RunningJobCounter = hogGroupConcurrencyCounters(hogGroupNumber)

      if (c.getMax == maxConcurrencyPerHogGroup) { exactlyAtLimit += 1 } else {
        Assertions.assert(c.getMax <= maxConcurrencyPerHogGroup, s"(asserting maxActualConcurrency for each hog group <= maxRequestedConcurrency per hog group)")
        Assertions.assert(c.getMax >= maxConcurrencyPerHogGroup * 0.95, s"(asserting maxActualConcurrency for each hog group >= (95% of maxRequestedConcurrency per hog group))")
      }
    }
    Assertions.assert(exactlyAtLimit >= (totalHogGroups * 0.95), "(at least 95% of the hog groups reached the full concurrency limit)")

    workflows foreach { system.stop(_) }
    system.stop(tokenDispenserUnderTest)
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
