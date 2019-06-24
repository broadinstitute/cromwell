package wom.util

import com.typesafe.config.ConfigException.BadValue
import com.typesafe.config.ConfigFactory
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.NonNegative
import eu.timepit.refined.refineMV
import io.circe.Json
import net.ceedubs.ficus.Ficus._
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{EitherValues, FlatSpec, Matchers}

class YamlUtilsSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks with EitherValues {

  behavior of "YamlUtils"

  private val illegalYamlTests = Table(
    ("description", "yaml", "maxNodes", "exceptionMessage"),
    (
      // https://stackoverflow.com/a/16986339
      "a stackoverflow yaml bomb",
      "{a: &b [*b]}",
      refineMV[NonNegative](10),
      "Loop detected",
    ),
    (
      // https://bitbucket.org/asomov/snakeyaml/wiki/Documentation#markdown-header-aliases
      "a recursive yaml example from the snakeyaml wiki",
      "&A [ *A ]",
      refineMV[NonNegative](10),
      "Loop detected",
    ),
    (
      "a recursive yaml based on mappings back to the root",
      """|a: &b
         |  c: *b
         |""".stripMargin,
      refineMV[NonNegative](10),
      "Loop detected",
    ),
    (
      // https://bitbucket.org/asomov/snakeyaml-engine/src/41b3845/src/test/resources/recursive/recursive-set-1.yaml
      "a recursive yaml set",
      """|# Sets are represented as a
         |# Mapping where each key is
         |# associated with a null value
         |--- !!set &key
         |? Mark McGwire
         |? Sammy Sosa
         |? *key
         |""".stripMargin,
      refineMV[NonNegative](10),
      "Loop detected",
    ),
    (
      // https://en.wikipedia.org/w/index.php?title=Billion_laughs_attack&oldid=871224525#Variations
      "a billion laughs when limited to 10,000 nodes",
      """|a: &a ["lol","lol","lol","lol","lol","lol","lol","lol","lol"]
         |b: &b [*a,*a,*a,*a,*a,*a,*a,*a,*a]
         |c: &c [*b,*b,*b,*b,*b,*b,*b,*b,*b]
         |d: &d [*c,*c,*c,*c,*c,*c,*c,*c,*c]
         |e: &e [*d,*d,*d,*d,*d,*d,*d,*d,*d]
         |f: &f [*e,*e,*e,*e,*e,*e,*e,*e,*e]
         |g: &g [*f,*f,*f,*f,*f,*f,*f,*f,*f]
         |h: &h [*g,*g,*g,*g,*g,*g,*g,*g,*g]
         |i: &i [*h,*h,*h,*h,*h,*h,*h,*h,*h]
         |""".stripMargin,
      refineMV[NonNegative](10000),
      "Loop detection halted at 10000 nodes",
    ),
    (
      "a null yaml",
      null,
      refineMV[NonNegative](0),
      null,
    ),
    (
      "an empty yaml",
      "",
      refineMV[NonNegative](1),
      "null",
    ),
    (
      "an empty yaml mapping when limited to zero nodes",
      "{}",
      refineMV[NonNegative](0),
      "Loop detection halted at 0 nodes",
    ),
    (
      "an empty yaml sequence when limited to zero nodes",
      "[]",
      refineMV[NonNegative](0),
      "Loop detection halted at 0 nodes",
    ),
    (
      "a yaml without a closing brace",
      "{",
      refineMV[NonNegative](10),
      """|while parsing a flow node
         | in 'reader', line 1, column 2:
         |    {
         |     ^
         |expected the node content, but found '<stream end>'
         | in 'reader', line 1, column 2:
         |    {
         |     ^
         |""".stripMargin,
    ),
  )

  private val legalYamlTests = Table(
    ("description", "yaml", "maxNodes", "expectedJson"),
    (
      "an empty yaml mapping when limited to one node",
      "{}",
      refineMV[NonNegative](1),
      Json.obj(),
    ),
    (
      "an empty yaml sequence when limited to one node",
      "[]",
      refineMV[NonNegative](1),
      Json.arr(),
    ),
    (
      "a yaml with the same node for a key and value",
      """|? &a b
         |: *a
         |""".stripMargin,
      refineMV[NonNegative](3),
      Json.obj("b" -> Json.fromString("b")),
    ),
    (
      "a yaml with the same nodes in a sequence",
      "[ &a b, *a ]",
      refineMV[NonNegative](3),
      Json.arr(Json.fromString("b"), Json.fromString("b")),
    ),
  )

  forAll(illegalYamlTests) { (description, yaml, maxNodes, exceptionMessage) =>
    it should s"fail to parse $description" in {
      val exception: Exception = YamlUtils.parse(yaml, maxNodes).left.value
      exception should have message exceptionMessage
    }
  }

  forAll(legalYamlTests) { (description, yaml, maxNodes, expectedJson) =>
    it should s"parse $description" in {
      YamlUtils.parse(yaml, maxNodes) should be(Right(expectedJson))
    }
  }

  it should "fail to parse a billion laughs with the default configuration" in {
    val yaml =
      """|a: &a ["lol","lol","lol","lol","lol","lol","lol","lol","lol"]
         |b: &b [*a,*a,*a,*a,*a,*a,*a,*a,*a]
         |c: &c [*b,*b,*b,*b,*b,*b,*b,*b,*b]
         |d: &d [*c,*c,*c,*c,*c,*c,*c,*c,*c]
         |e: &e [*d,*d,*d,*d,*d,*d,*d,*d,*d]
         |f: &f [*e,*e,*e,*e,*e,*e,*e,*e,*e]
         |g: &g [*f,*f,*f,*f,*f,*f,*f,*f,*f]
         |h: &h [*g,*g,*g,*g,*g,*g,*g,*g,*g]
         |i: &i [*h,*h,*h,*h,*h,*h,*h,*h,*h]
         |""".stripMargin
    val exception: Exception = YamlUtils.parse(yaml).left.value
    exception should have message "Loop detection halted at 1000000 nodes" // <-- Updating here? Also get the docs too!
  }

  it should "fail to parse deeply nested yaml sequence with the default configuration" in {
    val nesting = 1000000
    val yaml = ("[" * nesting) + ("]" * nesting)
    val exception: Exception = YamlUtils.parse(yaml).left.value
    exception should have message "Parsing halted at node depth 1000" // <-- Updating here? Also get the docs too!
  }

  it should "fail to parse deeply nested yaml mapping with the default configuration" in {
    val nesting = 1000000
    val yaml = ("{a:" * nesting) + "b" + ("}" * nesting)
    val exception: Exception = YamlUtils.parse(yaml).left.value
    exception should have message "Parsing halted at node depth 1000" // <-- Updating here? Also get the docs too!
  }

  it should "not parse a config with a negative value" in {
    import wom.util.YamlUtils.refinedNonNegativeReader
    the[BadValue] thrownBy {
      ConfigFactory.parseString("non-negative: -1").as[Int Refined NonNegative]("non-negative")
    } should have message "Invalid value at 'non-negative': Predicate (-1 < 0) did not fail."
  }
}
