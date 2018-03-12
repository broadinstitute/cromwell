package cromwell.engine.workflow.tokens

import akka.testkit.TestProbe
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.core.TestKitSuite
import org.scalatest.{FlatSpecLike, Matchers}

class RoundRobinQueueIteratorSpec extends TestKitSuite with FlatSpecLike with Matchers {
  behavior of "RoundRobinQueueIterator"
  
  val InfiniteTokenType = JobExecutionTokenType("infinite", None)
  val Pool1 = JobExecutionTokenType("pool1", Option(1))
  val Pool2 = JobExecutionTokenType("pool2", Option(2))
  
  it should "be empty if there's no queue" in {
    RoundRobinQueueIterator(List.empty).hasNext shouldBe false
  }

  it should "return an element if a queue can dequeue" in {
    val probe1 = TestProbe().ref
    val queues = List(
      TokenQueue(InfiniteTokenType).enqueue(probe1)
    )
    val iterator = RoundRobinQueueIterator(queues)
    iterator.hasNext shouldBe true
    val next = iterator.next()
    next.actor shouldBe probe1
    next.lease.get().jobExecutionTokenType shouldBe InfiniteTokenType
    iterator.hasNext shouldBe false

    iterator.updatedQueues.head.queue shouldBe empty 
  }
  
  it should "rotate between queues" in {
    val probe1 = TestProbe("probe-1").ref
    val probe2 = TestProbe("probe-2").ref
    val probe3 = TestProbe("probe-3").ref
    val queues = List(
      TokenQueue(InfiniteTokenType).enqueue(probe1).enqueue(probe3),
      TokenQueue(Pool2).enqueue(probe2)
    )
    val iterator = RoundRobinQueueIterator(queues)
    iterator.hasNext shouldBe true
    
    var next = iterator.next()
    next.actor shouldBe probe1
    next.lease.get().jobExecutionTokenType shouldBe InfiniteTokenType
    iterator.hasNext shouldBe true
    
    next = iterator.next()
    next.actor shouldBe probe2
    next.lease.get().jobExecutionTokenType shouldBe Pool2
    iterator.hasNext shouldBe true

    next = iterator.next()
    next.actor shouldBe probe3
    next.lease.get().jobExecutionTokenType shouldBe InfiniteTokenType
    iterator.hasNext shouldBe false

    iterator.updatedQueues.head.queue shouldBe empty
  }

  it should "respect pool sizes" in {
    val probe1 = TestProbe("probe-1").ref
    val probe2 = TestProbe("probe-2").ref
    val probe3 = TestProbe("probe-3").ref
    val probe4 = TestProbe("probe-4").ref
    val probe5 = TestProbe("probe-5").ref
    val queues = List(
      TokenQueue(Pool1)
        .enqueue(probe1)
        .enqueue(probe3),
      TokenQueue(Pool2)
        .enqueue(probe2)
        .enqueue(probe4)
        .enqueue(probe5)
    )
    val iterator = RoundRobinQueueIterator(queues)
    iterator.hasNext shouldBe true

    var next = iterator.next()
    next.actor shouldBe probe1
    val probe1Lease = next.lease
    probe1Lease.get().jobExecutionTokenType shouldBe Pool1
    iterator.hasNext shouldBe true

    next = iterator.next()
    next.actor shouldBe probe2
    next.lease.get().jobExecutionTokenType shouldBe Pool2
    iterator.hasNext shouldBe true
    
    // Pool 1 has given its only token, so take from pool2 again
    next = iterator.next()
    next.actor shouldBe probe4
    next.lease.get().jobExecutionTokenType shouldBe Pool2
    
    // both queues still have actors but the pools are empty
    iterator.hasNext shouldBe false

    iterator.updatedQueues.head.queue.toList should contain theSameElementsAs List(probe3)
    iterator.updatedQueues.last.queue.toList should contain theSameElementsAs List(probe5)

    // If we release the first lease, we should be able to iterate on more time
    probe1Lease.release()
    iterator.hasNext shouldBe true
    next = iterator.next()
    next.actor shouldBe probe3
    next.lease.get().jobExecutionTokenType shouldBe Pool1
    iterator.hasNext shouldBe false

    iterator.updatedQueues.head.queue shouldBe empty
    iterator.updatedQueues.last.queue.toList should contain theSameElementsAs List(probe5)
  }
  
  it should "throw an exception when calling next() on an empty iterator" in {
    assertThrows[IllegalStateException](RoundRobinQueueIterator(List.empty).next())
  }
}
