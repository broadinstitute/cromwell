package cromwell.util

//import io.circe._
import JsonEditor._
import cats.data.NonEmptyList
import io.circe.Json
import io.circe.parser._
import org.scalatest.{FlatSpec, Matchers}
import cats.syntax.either._
//import org.scalatest.LoneElement._

class CirceTest extends FlatSpec with Matchers{

  val rawJson =
    """
        { "foo": "bar", "other":"baz"}
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

  it should "keep includes" in {
    val either = testJsonAndGetKeys(includeExcludeJson(_, Some(NonEmptyList.one("foo")), None))
    assert(either.right.get.head === "foo")
  }

  "Json Munger" should "remove excludes" in {
    val either = testJsonAndGetKeys(includeExcludeJson(_, None, Some(NonEmptyList.one("foo"))))
    assert(either.right.get.head === "other")
  }
}
