package centaur.test.metadata

import centaur.test.Operations
import org.scalatest.{FlatSpec, Matchers}
import spray.json._

class ExtractJobManagerStyleMetadataFieldsSpec extends FlatSpec with Matchers {
  behavior of "extracting Job Manager style metadata fields"

  it should "preserve the right set of fields, including call caching hits, correctly" in {
    val originalMetadata =
      """{
        |	"id": "blah",
        |	"unwanted": "field",
        |	"calls": {
        |		"wf.foo": [{
        |			"attempt": 1,
        |			"shardIndex": -1,
        |			"unwanted": "field",
        |			"callCaching": {
        |				"unwanted": "field",
        |				"hit": "true",
        |				"hitFailures": [{
        |					"blah blah": "blah"
        |				}]
        |			}
        |		}]
        |	}
        |}""".stripMargin

    val expectedExpectedMetadata =
      """{
        |	"id": "blah",
        |	"calls": {
        |		"wf.foo": [{
        |			"attempt": 1,
        |			"shardIndex": -1,
        |			"callCaching": {
        |				"hit": "true"
        |			}
        |		}]
        |	}
        |}""".stripMargin

    Operations.extractJmStyleMetadataFields(originalMetadata.parseJson.asJsObject) should be(expectedExpectedMetadata.parseJson.asJsObject)
  }
}
