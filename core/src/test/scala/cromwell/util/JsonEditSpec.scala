package cromwell.util

import JsonEditor._
import cats.data.NonEmptyList
import io.circe.Json
import io.circe.parser._
import org.scalatest.{FlatSpec, Matchers}
import cats.syntax.either._

class JsonEditSpec extends FlatSpec with Matchers{

  val rawJson =
    """
      {
        "foo": "bar",
         "other":"baz",
         "nested": {
           "inner": {
             "deep":"some",
             "keepme": "more",
             "wildcard": "again"
           }
         }
       }
      """.stripMargin

  val jsonEither: Either[String, Json] = parse(rawJson).leftMap(_.toString)

  def testJson(f: Json => Json): Either[String, Json] =
    for {
      json <- jsonEither
      newJson = f(json)
    }  yield newJson

  def testJsonAndGetKeys(f: Json => Json): Either[String, Iterable[String]] = {
    for {
      newJson <- testJson(f)
      keys <- newJson.hcursor.keys.toRight("no keys found!")
    }  yield keys
  }


  "Json Munger" should "remove excludes" in {
      val either = testJsonAndGetKeys(includeExcludeJson(_, None, Some(NonEmptyList.one("foo"))))
      assert(either.right.get.head === "other")
    }

  it should "remove nested keys excludes" in {
    val either = testJson(excludeJson(_, NonEmptyList.one("deep")))
    assert(either.right.get.hcursor.downField("nested").downField("inner").keys.get.head === "keepme")
  }

  it should "remove multiple nested keys excludes" in {
    val either = testJson(excludeJson(_, NonEmptyList.of("deep", "wildcard")))
    val keys = either.right.get.hcursor.downField("nested").downField("inner").keys.get
    assert(keys.head === "keepme")
    assert(keys.size === 1)
  }

  it should "keep includes" in {
    val either = testJsonAndGetKeys(includeExcludeJson(_, Some(NonEmptyList.one("foo")), None))
    assert(either.right.get.head === "foo")
  }

  it should "keep nested includes" in {
    val either = testJson(includeJson(_, NonEmptyList.one("keepme")))
    assert(either.right.get.hcursor.downField("nested").downField("inner").keys.get.head === "keepme")
  }

  it should "keep multiple nested includes" in {
    val either = testJson(includeJson(_, NonEmptyList.of("keepme", "wildcard")))
    val keys = either.right.get.hcursor.downField("nested").downField("inner").keys.get
    assert(keys.head === "keepme")
    assert(keys.tail.head === "wildcard")
  }
}
