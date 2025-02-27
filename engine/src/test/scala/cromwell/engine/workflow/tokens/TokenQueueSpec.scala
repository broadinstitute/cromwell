package cromwell.engine.workflow.tokens

import akka.testkit.TestProbe
import cromwell.core.JobToken.JobTokenType
import cromwell.core.TestKitSuite
import cromwell.engine.workflow.tokens.TokenQueue.TokenQueuePlaceholder
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class TokenQueueSpec extends TestKitSuite with AnyFlatSpecLike with Matchers {
  behavior of "TokenQueue"
  val tokenType = JobTokenType("pool1", Option(1), 1)

  val tokenEventLogger = NullTokenEventLogger

  it should "enqueue" in {
    val probe = TestProbe().ref
    val tq = TokenQueue(tokenType, tokenEventLogger)
    tq.queues shouldBe empty
    val queueWith1 = tq.enqueue(TokenQueuePlaceholder(probe, "hogGroupA"))
    queueWith1.queues.flatMap(_._2).head.actor shouldBe probe
    queueWith1.queueOrder should be(List("hogGroupA"))
    queueWith1.size shouldBe 1
    queueWith1.tokenType shouldBe tokenType
  }

  it should "dequeue" in {
    val probe = TestProbe().ref
    val tq = TokenQueue(tokenType, tokenEventLogger).enqueue(TokenQueuePlaceholder(probe, "hogGroupA"))
    tq.queueOrder should be(List("hogGroupA"))
    val dequeued = tq.dequeue(List.empty)
    dequeued.leasedActor shouldBe defined
    dequeued.leasedActor.get.actor shouldBe probe
    dequeued.leasedActor.get.lease.get().jobTokenType shouldBe tokenType
    dequeued.tokenQueue.queues shouldBe empty
    dequeued.tokenQueue.queueOrder should be(List.empty)
    dequeued.tokenQueue.size shouldBe 0
  }

  it should "dequeue only from non-exhausted group" in {
    val probe1 = TestProbe().ref
    val probe2 = TestProbe().ref
    val probe3 = TestProbe().ref
    val probe4 = TestProbe().ref
    val quotaExhaustedGroup = "groot-hog-group"
    val nonQuotaExhaustedGroup = "not-groot-hog-group"
    val quotaExhaustedGroups = List(quotaExhaustedGroup)
    val mockTokenType = JobTokenType("pool1", Option(10), 3)

    val tq1 = TokenQueue(mockTokenType, tokenEventLogger)
      .enqueue(TokenQueuePlaceholder(probe1, quotaExhaustedGroup))
      .enqueue(TokenQueuePlaceholder(probe2, nonQuotaExhaustedGroup))

    // queue order should still be first-come-first-serve
    tq1.queueOrder should be(List(quotaExhaustedGroup, nonQuotaExhaustedGroup))

    // only request from non-exhausted quota should be dequeued
    tq1.available(quotaExhaustedGroups) shouldBe true
    val dequeued1 = tq1.dequeue(quotaExhaustedGroups)
    dequeued1.leasedActor shouldBe defined
    dequeued1.leasedActor.get.actor shouldBe probe2
    dequeued1.leasedActor.get.lease.get().jobTokenType shouldBe mockTokenType

    val updatedTQ1 = dequeued1.tokenQueue

    // there are no requests from non-exhausted group in queue
    updatedTQ1.available(quotaExhaustedGroups) shouldBe false

    // queue should still contain request from quota exhausted group
    updatedTQ1.queues.size shouldBe 1
    updatedTQ1.queueOrder should be(List(quotaExhaustedGroup))

    // another 2 requests from non-exhausted group
    val tq2 = updatedTQ1
      .enqueue(TokenQueuePlaceholder(probe3, nonQuotaExhaustedGroup))
      .enqueue(TokenQueuePlaceholder(probe4, nonQuotaExhaustedGroup))

    tq2.queueOrder should be(List(quotaExhaustedGroup, nonQuotaExhaustedGroup))

    // only requests from non-exhausted quota should be dequeued
    tq2.available(quotaExhaustedGroups) shouldBe true
    // it should dequeue probe3
    val dequeued2 = tq2.dequeue(quotaExhaustedGroups)
    dequeued2.leasedActor shouldBe defined
    dequeued2.leasedActor.get.actor shouldBe probe3
    dequeued2.leasedActor.get.lease.get().jobTokenType shouldBe mockTokenType

    val tq3 = dequeued2.tokenQueue

    // queue should contain requests from 2 groups for probe1 and probe4
    tq3.queues.size shouldBe 2
    tq3.queueOrder should be(List(quotaExhaustedGroup, nonQuotaExhaustedGroup))

    // 'groot-hog-group' is no longer in quota exhaustion

    // next dequeue should be from 'groot-hog-group' for probe1
    tq3.available(List.empty) shouldBe true
    // it should dequeue probe1
    val dequeued3 = tq3.dequeue(List.empty)
    dequeued3.leasedActor shouldBe defined
    dequeued3.leasedActor.get.actor shouldBe probe1
    dequeued3.leasedActor.get.lease.get().jobTokenType shouldBe mockTokenType

    val tq4 = dequeued3.tokenQueue

    // 1 request from probe4 should still exist and should be dequeued next
    tq4.available(List.empty) shouldBe true
    val dequeued4 = tq4.dequeue(List.empty)
    dequeued4.leasedActor shouldBe defined
    dequeued4.leasedActor.get.actor shouldBe probe4
    dequeued4.leasedActor.get.lease.get().jobTokenType shouldBe mockTokenType

    // all requests are fulfilled
    dequeued4.tokenQueue.queues shouldBe empty
  }

  it should "return true for available iff there's an element in the queue and room in the pool" in {
    val probe1 = TestProbe().ref
    val probe2 = TestProbe().ref
    val tq = TokenQueue(tokenType, tokenEventLogger)
    // queue is empty
    tq.available(List.empty) shouldBe false

    // with something in the queue, available is true
    val queueWith1 = tq.enqueue(TokenQueuePlaceholder(probe1, "hogGroupA"))
    queueWith1.available(List.empty) shouldBe true

    val queueWith2 = queueWith1.enqueue(TokenQueuePlaceholder(probe2, "hogGroupA"))

    val dequeued = queueWith2.dequeue(List.empty)
    // probe 2 is still in there
    dequeued.tokenQueue.size shouldBe 1
    // pool is empty though so available should be false
    dequeued.tokenQueue.available(List.empty) shouldBe false

    // If we release the token, we should be available again
    dequeued.leasedActor.get.lease.release()
    dequeued.tokenQueue.available(List.empty) shouldBe true
  }

  it should "respect hog groups" in {
    val refA1 = TestProbe("A1").ref

    // Start with an empty queue which can dish out up to 1 token each, to up to 2 hog groups:
    val tq = TokenQueue(JobTokenType("pool1", Option(2), 2), tokenEventLogger)
    tq.available(List.empty) shouldBe false

    val queueWith1 = tq.enqueue(TokenQueuePlaceholder(refA1, "hogGroupA"))
    queueWith1.queueOrder should be(List("hogGroupA"))
    queueWith1.available(List.empty) shouldBe true

    val dequeueResult = queueWith1.dequeue(List.empty)
    dequeueResult.leasedActor.map(_.actor) should be(Some(refA1))
    val tokenDispensedAlreadyQueue = dequeueResult.tokenQueue

    // Since it's been dispensed, there's now nothing available in this queue:
    tokenDispensedAlreadyQueue.queueOrder should be(List.empty)
    tokenDispensedAlreadyQueue.available(List.empty) shouldBe false

    // Adding another entry in hogGroupA should not help:
    val refA2 = TestProbe("A2").ref
    val queueWith2 = tokenDispensedAlreadyQueue.enqueue(TokenQueuePlaceholder(refA2, "hogGroupA"))
    queueWith2.available(List.empty) shouldBe false

    // An entry in hogGroupB however *should* register as available:
    val refB1 = TestProbe("B1").ref
    val queueWithB1 = queueWith2.enqueue(TokenQueuePlaceholder(refB1, "hogGroupB"))
    queueWithB1.available(List.empty) shouldBe true

    // And the 'B' actor should be the one to be dequeued next:
    val dequeueResult2 = queueWithB1.dequeue(List.empty)
    dequeueResult2.leasedActor.map(_.actor) should be(Some(refB1))
    val tokenB1DispensedAlreadyQueue = dequeueResult2.tokenQueue

    // And again, now we don't have anything available:
    tokenB1DispensedAlreadyQueue.available(List.empty) shouldBe false

    // Adding another entry in hogGroupB should not help:
    val refB2 = TestProbe("B2").ref
    val queueWithB2 = tokenB1DispensedAlreadyQueue.enqueue(TokenQueuePlaceholder(refB2, "hogGroupB"))
    queueWithB2.available(List.empty) shouldBe false

    // Because the pool is now full, adding a third hog group won't make anything available either:
    val refC1 = TestProbe("C1").ref
    val queueWithA2B2andC1 = queueWithB2.enqueue(TokenQueuePlaceholder(refC1, "hogGroupC"))
    queueWithA2B2andC1.available(List.empty) shouldBe false

    // When we release B1's lease, we get availability again:
    dequeueResult2.leasedActor.foreach(_.lease.release())
    queueWithA2B2andC1.available(List.empty) shouldBe true

    // The next entry to be dequeued should be B2, not C1:
    val dequeueResult3 = queueWithA2B2andC1.dequeue(List.empty)
    dequeueResult3.leasedActor.map(_.actor) should be(Some(refB2))
    val queueWithA2andC1 = dequeueResult3.tokenQueue
    queueWithA2andC1.available(List.empty) should be(false)

    // But when B2 finishes, we do indeed get C1:
    dequeueResult3.leasedActor foreach { _.lease.release() }
    queueWithA2andC1.available(List.empty) should be(true)
    val dequeueResult4 = queueWithA2andC1.dequeue(List.empty)
    dequeueResult4.leasedActor.map(_.actor) should be(Some(refC1))
  }

  it should "enqueue and dequeue with multiple hog groups" in {
    import cromwell.engine.workflow.tokens.large.TokenDispenserBenchmark.{fillQueue, useEntireAvailability}

    val poolSize = 52
    val jobCount = 150
    val jobsAtATime = 5
    val hogFactor = 10
    val tokenType = JobTokenType("BIG_PAPI", Some(poolSize), hogFactor)

    val hogGroupCount = 25
    val hogGroups = (0 until hogGroupCount) map { i => s"hogGroup$i" }

    val jobsPerHogGroup = jobCount / hogGroupCount

    val fullQueue = fillQueue(TokenQueue(tokenType, tokenEventLogger), jobsPerHogGroup, hogGroups.toList)
    val actualFullQueueSizes = fullQueue.queueOrder.map(queueName => fullQueue.queues(queueName).size)
    val expectedFullQueueSizes = (0 until 25).toVector.map(_ => 6)
    actualFullQueueSizes should be(expectedFullQueueSizes)
    fullQueue.size should be(jobCount)
    fullQueue.queueOrder should be(hogGroups.toVector)

    val usedQueue = useEntireAvailability(fullQueue, jobsAtATime)
    val tokenQueueSizes = usedQueue.queueOrder.map(queueName => usedQueue.queues(queueName).size)

    // Most hog group queues have 4 jobs left, but the final two have 3 (because there are 25 hog groups and 52 tokens)
    val expectedQueueSizes = (0 until 23).toVector.map(_ => 4) ++ Vector(3, 3)
    tokenQueueSizes should be(expectedQueueSizes)
    // Make sure the hog group queues are cycling correctly (0 and 1 are at the end because they already had the extra tokens):
    val expectedOrder = (0 until 23).toVector.map(i => s"hogGroup${i + 2}") ++ Vector("hogGroup0", "hogGroup1")
    usedQueue.queueOrder should be(expectedOrder)

    usedQueue.size should be(jobCount - poolSize)

  }

}
