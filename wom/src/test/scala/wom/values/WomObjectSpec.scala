package wom.values

import org.scalatest.TryValues
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wom.types._

object WomObjectSpec {
  implicit class Trimmer(val string: String) extends AnyVal {
    // Return both the original version of the string as well as a version trimmed of a newline.
    def withTrimmed: List[String] = List(string, string.trim)
  }
}

class WomObjectSpec extends AnyFlatSpec with Matchers with TryValues {
  import WomObjectSpec._
  val correctTSV = "one\ttwo\tthree\tfour\none\tfour\tnine\tsixteen\n"
  val emptyTSV = ""
  val oneRowTSV = "one\ttwo\tthree\tfour\n"
  val nonHomogeneousTSV = "one\ttwo\tthree\none\ttwo\n"
  val arrayTSV = correctTSV + "one\teight\ttwentyseven\tsixtyfour\n"

  it should "read an Object from a correct TSV file" in {
    // Test both a version of the TSV with and without a trailing newline.
    correctTSV.withTrimmed foreach { tsv =>
      val parsed = WomObject.fromTsv(tsv)
      parsed should be a 'success
      val array: Array[WomObject] = parsed.success.value
      array should have size 1

      //Attributes
      array.head.values should contain key "one"
      array.head.values should contain key "two"
      array.head.values should contain key "three"
      array.head.values should contain key "four"

      //Values
      array.head.values.get("one") shouldBe Some(WomString("one"))
      array.head.values.get("two") shouldBe Some(WomString("four"))
      array.head.values.get("three") shouldBe Some(WomString("nine"))
      array.head.values.get("four") shouldBe Some(WomString("sixteen"))
    }
  }

  it should "NOT read from a TSV file with less than 2 rows" in {
    for {
      tsv <- List(emptyTSV, oneRowTSV)
      t <- tsv.withTrimmed
      _ = WomObject.fromTsv(t) should be a 'failure
    } yield ()
  }

  it should "NOT read from a non homogeneous TSV file" in {
    nonHomogeneousTSV.withTrimmed foreach {
      WomObject.fromTsv(_) should be a 'failure
    }
  }

  it should "serialize to TSV" in {
    correctTSV.withTrimmed foreach { tsv =>
      val obj = WomObject.fromTsv(tsv).get.head
      val serialized = obj.tsvSerialize
      serialized should be a 'success
      serialized.success.value shouldEqual correctTSV
    }
  }

  it should "read a WomArray[WomObject] from a correct TSV file" in {
    List(arrayTSV, arrayTSV.trim) foreach { tsv =>
      val parsed = WomObject.fromTsv(tsv)
      parsed should be a 'success
      val array: Array[WomObject] = parsed.success.value
      array should have size 2

      //Attributes
      array foreach { a =>
        a.values should contain key "one"
        a.values should contain key "two"
        a.values should contain key "three"
        a.values should contain key "four"
      }

      //Values
      array.head.values.get("one") shouldBe Some(WomString("one"))
      array.head.values.get("two") shouldBe Some(WomString("four"))
      array.head.values.get("three") shouldBe Some(WomString("nine"))
      array.head.values.get("four") shouldBe Some(WomString("sixteen"))

      array(1).values.get("one") shouldBe Some(WomString("one"))
      array(1).values.get("two") shouldBe Some(WomString("eight"))
      array(1).values.get("three") shouldBe Some(WomString("twentyseven"))
      array(1).values.get("four") shouldBe Some(WomString("sixtyfour"))
    }
  }

  it should "serialize a WomArray[WomObject] to TSV" in {
    List(arrayTSV, arrayTSV.trim) foreach { tsv =>
      val array = WomArray(WomArrayType(WomObjectType), WomObject.fromTsv(tsv).get)
      val serialized = array.tsvSerialize
      serialized should be a 'success
      serialized.success.value shouldEqual arrayTSV
    }
  }
}
