package cromwell.engine.workflow.tokens

import akka.testkit.TestProbe
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.core.TestKitSuite
import org.scalatest.{FlatSpecLike, Matchers}

class TokenQueueSpec extends TestKitSuite with FlatSpecLike with Matchers {
  behavior of "TokenQueue"
  val tokenType = JobExecutionTokenType("pool1", Option(1))
  
  it should "enqueue" in {
    val probe = TestProbe().ref
    val tq = TokenQueue(tokenType)
    tq.queue shouldBe empty
    val queueWith1 = tq.enqueue(probe)
    queueWith1.queue.head shouldBe probe
    queueWith1.size shouldBe 1
    queueWith1.tokenType shouldBe tokenType
  }

  it should "deenqueue" in {
    val probe = TestProbe().ref
    val tq = TokenQueue(tokenType)
      .enqueue(probe)
    val dequeued = tq.dequeueOption
    dequeued shouldBe defined
    dequeued.get.leasedActor.actor shouldBe probe
    dequeued.get.leasedActor.lease.get().jobExecutionTokenType shouldBe tokenType
    dequeued.get.tokenQueue.queue shouldBe empty
    dequeued.get.tokenQueue.size shouldBe 0
  }
  
  it should "return true for available iff there's an element in the queue and room in the pool" in {
    val probe1 = TestProbe().ref
    val probe2 = TestProbe().ref
    val tq = TokenQueue(tokenType)
    // queue is empty
    tq.available shouldBe false
    
    // with something in the queue, available is true
    val queueWith1 = tq.enqueue(probe1)
    queueWith1.available shouldBe true

    val queueWith2 = queueWith1.enqueue(probe2)
    
    val dequeued = queueWith2.dequeueOption.get
    // probe 2 is still in there
    dequeued.tokenQueue.size shouldBe 1
    // pool is empty though so available should be false
    dequeued.tokenQueue.available shouldBe false

    // If we release the token, we should be available again
    dequeued.leasedActor.lease.release()
    dequeued.tokenQueue.available shouldBe true
  }
}
