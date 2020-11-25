package common.collections

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class WeightedQueueSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "WeightedQueue"

  it should "enqueue and dequeue elements" in {
    // A queue of strings for which the weight is the number of char in the string
    val q = WeightedQueue.empty[String, Int](_.length)
    val q2 = q.enqueue("hello")
    q2.weight shouldBe 5
    q2.innerQueue.size shouldBe 1
    val q3 = q2.enqueue("bye")
    q3.weight shouldBe 8
    q3.innerQueue.size shouldBe 2
    val (head, q4) = q3.dequeue
    head shouldBe "hello"
    q4.weight shouldBe 3
    q4.innerQueue.size shouldBe 1
  }

  it should "dequeue optionally" in {
    // A queue of strings for which the weight is the number of char in the string
    val q = WeightedQueue.empty[String, Int](_.length)
    val q2 = q.enqueue("hello")
    val dequeued1 = q2.dequeueOption
    dequeued1 shouldBe defined
    val (head, q3) = dequeued1.get
    head shouldBe "hello"
    q3.innerQueue shouldBe empty
    q3.dequeueOption shouldBe empty
  }

  it should "behead the queue" in {
    // A queue of strings for which the weight is the number of char in the string
    val q = WeightedQueue.empty[String, Int](_.length)
    val q2 = q.enqueue("hello")
      .enqueue("hola")
      .enqueue("bonjour")
    val (head, q3) = q2.behead(10)
    head shouldBe Vector("hello", "hola")
    q3.weight shouldBe 7
    q3.innerQueue.size shouldBe 1
  }

}
