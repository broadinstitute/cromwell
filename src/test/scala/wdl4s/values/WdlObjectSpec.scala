package wdl4s.values

import wdl4s.types.{WdlArrayType, WdlObjectType}
import org.scalatest.{FlatSpec, Matchers, TryValues}

object WdlObjectSpec {
  implicit class Trimmer(val string: String) extends AnyVal {
    // Return both the original version of the string as well as a version trimmed of a newline.
    def withTrimmed: List[String] = List(string, string.trim)
  }
}

class WdlObjectSpec extends FlatSpec with Matchers with TryValues {
  import WdlObjectSpec._
  val correctTSV = "one\ttwo\tthree\tfour\none\tfour\tnine\tsixteen\n"
  val emptyTSV = ""
  val oneRowTSV = "one\ttwo\tthree\tfour\n"
  val nonHomogeneousTSV = "one\ttwo\tthree\none\ttwo\n"
  val arrayTSV = correctTSV + "one\teight\ttwentyseven\tsixtyfour\n"

  it should "read an Object from a correct TSV file" in {
    // Test both a version of the TSV with and without a trailing newline.
    correctTSV.withTrimmed foreach { tsv =>
      val parsed = WdlObject.fromTsv(tsv)
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
  }

  it should "NOT read from a TSV file with less than 2 rows" in {
    for {
      tsv <- List(emptyTSV, oneRowTSV)
      t <- tsv.withTrimmed
      _ = WdlObject.fromTsv(t) should be a 'failure
    } yield ()
  }

  it should "NOT read from a non homogeneous TSV file" in {
    nonHomogeneousTSV.withTrimmed foreach {
      WdlObject.fromTsv(_) should be a 'failure
    }
  }

  it should "serialize to TSV" in {
    correctTSV.withTrimmed foreach { tsv =>
      val obj = WdlObject.fromTsv(tsv).get.head
      val serialized = obj.tsvSerialize
      serialized should be a 'success
      serialized.success.value shouldEqual correctTSV
    }
  }

  it should "read a WdlArray[WdlObject] from a correct TSV file" in {
    List(arrayTSV, arrayTSV.trim) foreach { tsv =>
      val parsed = WdlObject.fromTsv(tsv)
      parsed should be a 'success
      val array: Array[WdlObject] = parsed.success.value
      array should have size 2

      //Attributes
      array foreach { a =>
        a.value should contain key "one"
        a.value should contain key "two"
        a.value should contain key "three"
        a.value should contain key "four"
      }

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
  }

  it should "serialize a WdlArray[WdlObject] to TSV" in {
    List(arrayTSV, arrayTSV.trim) foreach { tsv =>
      val array = WdlArray(WdlArrayType(WdlObjectType), WdlObject.fromTsv(tsv).get)
      val serialized = array.tsvSerialize
      serialized should be a 'success
      serialized.success.value shouldEqual arrayTSV
    }
  }
}
