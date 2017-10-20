package lenthall.collections

import org.scalatest.{FlatSpec, Matchers}
import lenthall.collections.EnhancedCollections._

class EnhancedCollectionsSpec extends FlatSpec with Matchers {
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
}
