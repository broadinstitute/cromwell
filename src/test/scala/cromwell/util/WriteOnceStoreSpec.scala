package cromwell.util

import org.scalatest.{FlatSpec, Matchers}

class WriteOnceStoreSpec extends FlatSpec with Matchers {

  "A WriteOnceStore" should "allow you to insert a value" in {
    val writeOnceStore = new WriteOnceStore[String, Int]
    val output = writeOnceStore.insert("hi", 5)
    output.get.get("hi").get shouldEqual 5
  }

  it should "not allow you to insert the same value twice" in {
    val writeOnceStore = new WriteOnceStore[String, Int]
    writeOnceStore.insert("hi", 5)
    val blowUp = writeOnceStore.insert("hi", 5)
    blowUp should be a 'failure
  }

  it should "allow you to get a Map" in {
    val writeOnceStore = new WriteOnceStore[String, Int]
    writeOnceStore.insert("hi", 5)
    val theMap = writeOnceStore.toMap
    theMap shouldBe a [Map[_, _]]
    theMap.get("hi").get shouldEqual 5
  }
}
