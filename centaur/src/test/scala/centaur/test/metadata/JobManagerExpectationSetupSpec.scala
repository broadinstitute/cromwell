package centaur.test.metadata

import centaur.test.Operations
import org.scalatest.{FlatSpec, Matchers}
import spray.json._

class JobManagerExpectationSetupSpec extends FlatSpec with Matchers {
  behavior of "Job Manager expectation setup"

  it should "preserve fields, including call caching hits, correctly" in {
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

    Operations.setUpJmStyleMetadataExpectation(originalMetadata.parseJson.asJsObject) should be(expectedExpectedMetadata.parseJson.asJsObject)
  }
}
