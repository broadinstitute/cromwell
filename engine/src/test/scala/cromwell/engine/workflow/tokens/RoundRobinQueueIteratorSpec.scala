package cromwell.engine.workflow.tokens

import akka.testkit.TestProbe
import cromwell.core.JobToken.JobTokenType
import cromwell.core.TestKitSuite
import cromwell.engine.workflow.tokens.TokenQueue.TokenQueuePlaceholder
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class RoundRobinQueueIteratorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers {
  behavior of "RoundRobinQueueIterator"

  val InfiniteTokenType = JobTokenType("infinite", None, 1)
  val Pool1 = JobTokenType("pool1", Option(1), 1)
  val Pool2 = JobTokenType("pool2", Option(2), 1)

  val tokenEventLogger = NullTokenEventLogger

  it should "be empty if there's no queue" in {
    new RoundRobinQueueIterator(List.empty, 0, List.empty).hasNext shouldBe false
  }

  it should "return empty if all token requests are from quota exhausted groups" in {
    val probe1 = TestProbe().ref
    val probe2 = TestProbe().ref
    val probe3 = TestProbe().ref

    val quotaExhaustedGroup = "groot-hog-group"

    val queues = List(
      TokenQueue(InfiniteTokenType, tokenEventLogger)
        .enqueue(TokenQueuePlaceholder(probe1, quotaExhaustedGroup))
        .enqueue(TokenQueuePlaceholder(probe2, quotaExhaustedGroup))
        .enqueue(TokenQueuePlaceholder(probe3, quotaExhaustedGroup))
    )
    val iterator = new RoundRobinQueueIterator(queues, 0, List(quotaExhaustedGroup))

    iterator.hasNext shouldBe false
  }

  it should "return an element if a queue can dequeue (no quota exhausted groups)" in {
    val probe1 = TestProbe().ref
    val queues = List(
      TokenQueue(InfiniteTokenType, tokenEventLogger).enqueue(TokenQueuePlaceholder(probe1, "hogGroupA"))
    )
    val iterator = new RoundRobinQueueIterator(queues, 0, List.empty)
    iterator.hasNext shouldBe true
    val next = iterator.next()
    next.actor shouldBe probe1
    next.lease.get().jobTokenType shouldBe InfiniteTokenType
    iterator.hasNext shouldBe false

    iterator.updatedQueues.head.queues shouldBe empty
  }

  it should "return an element from queue whose hog group is not in quota exhausted groups" in {
    val probe1 = TestProbe().ref
    val probe2 = TestProbe().ref
    val probe3 = TestProbe().ref

    val quotaExhaustedGroup = "groot-hog-group"
    val nonQuotaExhaustedGroup = "rocket-hog-group"

    val queues = List(
      TokenQueue(InfiniteTokenType, tokenEventLogger)
        .enqueue(TokenQueuePlaceholder(probe1, quotaExhaustedGroup))
        .enqueue(TokenQueuePlaceholder(probe2, quotaExhaustedGroup))
        .enqueue(TokenQueuePlaceholder(probe3, nonQuotaExhaustedGroup))
    )
    val iterator = new RoundRobinQueueIterator(queues, 0, List(quotaExhaustedGroup))

    // there should be 1 element that can be dequeued
    iterator.hasNext shouldBe true

    val next = iterator.next()
    next.actor shouldBe probe3
    next.lease.get().jobTokenType shouldBe InfiniteTokenType

    // since 'groot-hog-group' is in quota exhausted state there is no more element to dequeue
    iterator.hasNext shouldBe false

    // the queue should not be empty as it still consists requests from quota exhausted hog group
    iterator.updatedQueues.head.queues.size shouldBe 1
    iterator.updatedQueues.head.queues.head._2.size shouldBe 2
  }

  it should "rotate between queues" in {
    val probe1 = TestProbe("probe-1").ref
    val probe2 = TestProbe("probe-2").ref
    val probe3 = TestProbe("probe-3").ref
    val queues = List(
      TokenQueue(InfiniteTokenType, tokenEventLogger)
        .enqueue(TokenQueuePlaceholder(probe1, "hogGroupA"))
        .enqueue(TokenQueuePlaceholder(probe3, "hogGroupA")),
      TokenQueue(Pool2, tokenEventLogger).enqueue(TokenQueuePlaceholder(probe2, "hogGroupA"))
    )
    val iterator = new RoundRobinQueueIterator(queues, 0, List.empty)
    iterator.hasNext shouldBe true

    var next = iterator.next()
    next.actor shouldBe probe1
    next.lease.get().jobTokenType shouldBe InfiniteTokenType
    iterator.hasNext shouldBe true

    next = iterator.next()
    next.actor shouldBe probe2
    next.lease.get().jobTokenType shouldBe Pool2
    iterator.hasNext shouldBe true

    next = iterator.next()
    next.actor shouldBe probe3
    next.lease.get().jobTokenType shouldBe InfiniteTokenType
    iterator.hasNext shouldBe false

    iterator.updatedQueues.head.queues shouldBe empty
  }

  it should "respect pool sizes" in {
    val probe1 = TestProbe("probe-1").ref
    val probe2 = TestProbe("probe-2").ref
    val probe3 = TestProbe("probe-3").ref
    val probe4 = TestProbe("probe-4").ref
    val probe5 = TestProbe("probe-5").ref
    val queues = List(
      TokenQueue(Pool1, tokenEventLogger)
        .enqueue(TokenQueuePlaceholder(probe1, "hogGroupA"))
        .enqueue(TokenQueuePlaceholder(probe3, "hogGroupA")),
      TokenQueue(Pool2, tokenEventLogger)
        .enqueue(TokenQueuePlaceholder(probe2, "hogGroupA"))
        .enqueue(TokenQueuePlaceholder(probe4, "hogGroupA"))
        .enqueue(TokenQueuePlaceholder(probe5, "hogGroupA"))
    )
    val iterator = new RoundRobinQueueIterator(queues, 0, List.empty)
    iterator.hasNext shouldBe true

    var next = iterator.next()
    next.actor shouldBe probe1
    val probe1Lease = next.lease
    probe1Lease.get().jobTokenType shouldBe Pool1
    iterator.hasNext shouldBe true

    next = iterator.next()
    next.actor shouldBe probe2
    next.lease.get().jobTokenType shouldBe Pool2
    iterator.hasNext shouldBe true

    // Pool 1 has given its only token, so take from pool2 again
    next = iterator.next()
    next.actor shouldBe probe4
    next.lease.get().jobTokenType shouldBe Pool2

    // both queues still have actors but the pools are empty
    iterator.hasNext shouldBe false

    iterator.updatedQueues.head.queues.flatMap(_._2).toList.map(_.actor) should contain theSameElementsAs List(probe3)
    iterator.updatedQueues.last.queues.flatMap(_._2).toList.map(_.actor) should contain theSameElementsAs List(probe5)

    // If we release the first lease, we should be able to iterate one more time
    probe1Lease.release()
    iterator.hasNext shouldBe true
    next = iterator.next()
    next.actor shouldBe probe3
    next.lease.get().jobTokenType shouldBe Pool1
    iterator.hasNext shouldBe false

    iterator.updatedQueues.head.queues.flatMap(_._2) shouldBe empty
    iterator.updatedQueues.last.queues.flatMap(_._2).toList.map(_.actor) should contain theSameElementsAs List(probe5)
  }

  it should "respect pool sizes and quota exhausted groups" in {
    val probe1 = TestProbe("probe-1").ref
    val probe2 = TestProbe("probe-2").ref
    val probe3 = TestProbe("probe-3").ref
    val probe4 = TestProbe("probe-4").ref
    val probe5 = TestProbe("probe-5").ref

    val quotaExhaustedGroup = "groot-hog-group"
    val nonQuotaExhaustedGroup = "rocket-hog-group"

    val queues = List(
      TokenQueue(Pool1, tokenEventLogger)
        .enqueue(TokenQueuePlaceholder(probe1, nonQuotaExhaustedGroup))
        .enqueue(TokenQueuePlaceholder(probe3, nonQuotaExhaustedGroup)),
      TokenQueue(InfiniteTokenType, tokenEventLogger)
        .enqueue(TokenQueuePlaceholder(probe2, quotaExhaustedGroup))
        .enqueue(TokenQueuePlaceholder(probe4, quotaExhaustedGroup))
        .enqueue(TokenQueuePlaceholder(probe5, quotaExhaustedGroup))
    )

    val iterator = new RoundRobinQueueIterator(queues, 0, List(quotaExhaustedGroup))

    // there are 2 requests for non-exhausted group and Pool1 has 1 token to lease
    iterator.hasNext shouldBe true

    // token should be leased to probe1 from non-exhausted group
    var next = iterator.next()
    next.actor shouldBe probe1
    val probe1Lease = next.lease
    probe1Lease.get().jobTokenType shouldBe Pool1

    // Pool1 has exhausted all its token and the requests from InfiniteTokenType all belong to quota exhausted group
    iterator.hasNext shouldBe false

    // if the first lease is released, we should be able to assign again from Pool1
    probe1Lease.release()
    iterator.hasNext shouldBe true
    next = iterator.next()
    next.actor shouldBe probe3
    next.lease.get().jobTokenType shouldBe Pool1

    // the only requests remain from quota exhausted group hence we can't iterate
    iterator.hasNext shouldBe false
  }

  it should "throw an exception when calling next() on an empty iterator" in {
    assertThrows[IllegalStateException](new RoundRobinQueueIterator(List.empty, 0, List.empty).next())
  }
}
