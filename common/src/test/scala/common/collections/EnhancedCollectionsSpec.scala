package common.collections

import common.collections.EnhancedCollections._
import org.scalatest.flatspec.AsyncFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.enablers.Emptiness._

import scala.collection.immutable.Queue

class EnhancedCollectionsSpec extends AsyncFlatSpec with Matchers {
  behavior of "EnhancedCollections"

  it should "filter a List by type and return a List" in {
    val objectList = List("hello", 3, None, "world")
    val stringList = objectList.filterByType[String]

    stringList should be(List("hello", "world"))
  }

  it should "filter a Set by type and return a Set" in {
    val objectSet = Set("hello", 3, None, "world")
    val intSet: Set[Int] = objectSet.filterByType[Int]

    intSet should be(Set(3))
  }

  it should "find the first Int in a List" in {
    val objectSet = List("hello", 3, None, 4, "world")
    objectSet.firstByType[Int] should be(Some(3))
  }

  behavior of "EnhancedQueue"
  
  it should "behead an empty queue" in {
    val q = Queue.empty[Int]
    val DeQueued(head, newQueue) = q.takeWhileWeighted(10, identity[Int], None)
    head shouldBe empty
    newQueue shouldBe empty
  }

  it should "behead a non empty queue with a max length of 0" in {
    val q = Queue(1)
    val DeQueued(head, newQueue) = q.takeWhileWeighted(10, identity[Int], Option(0))
    head shouldBe empty
    newQueue shouldBe q
  }

  it should "use the weight function to compute weight and respect the max weight" in {
    val q = Queue(1, 2, 5)
    val DeQueued(head, newQueue) = q.takeWhileWeighted(maxWeight = 10, _ * 2, None)
    head shouldBe Vector(1, 2) // weight is 5 (1 * 2 + 2 * 2)
    newQueue shouldBe Queue(5)
  }

  it should "respect the max length" in {
    val q = Queue(1, 2, 5)
    val DeQueued(head, newQueue) = q.takeWhileWeighted(maxWeight = 10, identity[Int], Option(2))
    head shouldBe Vector(1, 2)
    newQueue shouldBe Queue(5)
  }

  it should "not enable the strict flag by default and stop at the last element that goes over maxWeight" in {
    val q = Queue(1, 2, 5)
    val DeQueued(head, newQueue) = q.takeWhileWeighted(maxWeight = 4, identity[Int], Option(3))
    head shouldBe Vector(1, 2)
    newQueue shouldBe Queue(5)
  }

  it should "respect when the strict flag is enabled and drop messages over maxWeight" in {
    val q = Queue(1, 5, 2)
    val DeQueued(head, newQueue) = q.takeWhileWeighted(maxWeight = 4, identity[Int], Option(3), strict = true)
    head shouldBe Vector(1, 2)
    newQueue shouldBe empty
  }

  it should "return the first element as the head even if it's over max weight and strict is disabled" in {
    val q = Queue(5, 1, 2)
    val DeQueued(head, newQueue) = q.takeWhileWeighted(maxWeight = 4, identity[Int], Option(3))
    head shouldBe Vector(5)
    newQueue shouldBe Queue(1, 2)
  }

  it should "drop the first element if it's over max weight and strict is enabled" in {
    val q = Queue(5, 1, 2)
    val DeQueued(head, newQueue) = q.takeWhileWeighted(maxWeight = 4, identity[Int], Option(3), strict = true)
    head shouldBe Vector(1, 2)
    newQueue shouldBe empty
  }
}
