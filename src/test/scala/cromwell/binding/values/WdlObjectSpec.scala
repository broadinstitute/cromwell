package cromwell.binding.values

import cromwell.binding.types.{WdlArrayType, WdlObjectType}
import cromwell.engine.Hashing._
import cromwell.engine.backend.local.SharedFileSystem
import cromwell.util.HashUtil
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers, TryValues}

class WdlObjectSpec extends FlatSpec with Matchers with TryValues with HashUtil {

  val correctTSV = "one\ttwo\tthree\tfour\none\tfour\tnine\tsixteen"
  val emptyTSV = ""
  val oneRowTSV = "one\ttwo\tthree\tfour"
  val nonHomogeneousTS = "onet\ttwo\tthree\none\ttwo"
  val arrayTSV = correctTSV + "\none\teight\ttwentyseven\tsixtyfour"

  it should "read an Object from a correct TSV file" in {
    val parsed = WdlObject.fromTsv(correctTSV)
    parsed should be a 'success
    val array: Array[WdlObject] = parsed.success.value
    array should have size 1

    //Attributes
    array.head.value should contain key "one"
    array.head.value should contain key "two"
    array.head.value should contain key "three"
    array.head.value should contain key "four"

    //Values
    array.head.value.get("one") shouldBe Some(WdlString("one"))
    array.head.value.get("two") shouldBe Some(WdlString("four"))
    array.head.value.get("three") shouldBe Some(WdlString("nine"))
    array.head.value.get("four") shouldBe Some(WdlString("sixteen"))
  }

  it should "NOT read from a TSV file with less than 2 rows" in {
    WdlObject.fromTsv(emptyTSV) should be a 'failure
    WdlObject.fromTsv(oneRowTSV) should be a 'failure
  }

  it should "NOT read from a non homogeneous TSV file" in {
    WdlObject.fromTsv(nonHomogeneousTS) should be a 'failure
  }

  it should "serialize to TSV" in {
    val obj = WdlObject.fromTsv(correctTSV).get.head
    val serialized = obj.tsvSerialize
    serialized should be a 'success
    serialized.success.value shouldEqual correctTSV
  }

  it should "read a WdlArray[WdlObject] from a correct TSV file" in {
    val parsed = WdlObject.fromTsv(arrayTSV)
    parsed should be a 'success
    val array: Array[WdlObject] = parsed.success.value
    array should have size 2

    //Attributes
    array foreach { _.value should contain key "one" }
    array foreach { _.value should contain key "two" }
    array foreach { _.value should contain key "three" }
    array foreach { _.value should contain key "four" }

    //Values
    array.head.value.get("one") shouldBe Some(WdlString("one"))
    array.head.value.get("two") shouldBe Some(WdlString("four"))
    array.head.value.get("three") shouldBe Some(WdlString("nine"))
    array.head.value.get("four") shouldBe Some(WdlString("sixteen"))

    array(1).value.get("one") shouldBe Some(WdlString("one"))
    array(1).value.get("two") shouldBe Some(WdlString("eight"))
    array(1).value.get("three") shouldBe Some(WdlString("twentyseven"))
    array(1).value.get("four") shouldBe Some(WdlString("sixtyfour"))
  }

  it should "serialize a WdlArray[WdlObject] to TSV" in {
    val array = WdlArray(WdlArrayType(WdlObjectType), WdlObject.fromTsv(arrayTSV).get)
    val serialized = array.tsvSerialize
    serialized should be a 'success
    serialized.success.value shouldEqual arrayTSV
  }

  it should "produce correct Hash" in {
    implicit val fileHasher = SharedFileSystem.sharedFsFileHasher

    // Note that we don't need to test with every single possible WdlType as an attribute, because each type is responsible for producing a valid hash.

    val refObject = WdlObject(Map("key0" -> file1, "key1" -> string1))

    val emptyWdlObject: WdlObject = WdlObject(Map.empty)

    // Object with same hashes
    val sameHashTable = Table(
      ("obj1", "obj2"),

      // Just in case...
      (refObject, refObject),
      // Empty objects
      (emptyWdlObject, emptyWdlObject),
      // Same inner values but different WdlValues
      (refObject, WdlObject(Map("key0" -> sameAsfile1, "key1" -> sameAsString1))),
      // Shuffled ordering
      (refObject, WdlObject(Map("key1" -> sameAsString1, "key0" -> sameAsfile1)))
    )

    val differentHashTable = Table(
      ("obj1", "obj2"),

      (refObject, emptyWdlObject),
      // Missing an entry
      (refObject, WdlObject(Map("key0" -> file1))),
      // One more entry
      (refObject, WdlObject(Map("key0" -> file1, "key1" -> string1, "key2" -> string1))),
      // With a different key
      (refObject, WdlObject(Map("differentKey0" -> file1, "key1" -> string1))),
      // With a different value
      (refObject, WdlObject(Map("key0" -> anotherFile, "key1" -> string1)))
    )

    forAll(sameHashTable) { (obj1, obj2) => obj1.getHash should be(obj2.getHash) }
    forAll(differentHashTable) { (obj1, obj2) => obj1.getHash shouldNot be(obj2.getHash) }
  }

}
