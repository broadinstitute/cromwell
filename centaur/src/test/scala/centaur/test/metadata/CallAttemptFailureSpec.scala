package centaur.test.metadata

import centaur.test.metadata.CallAttemptFailureSpec._
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.JsonParser.ParsingException

class CallAttemptFailureSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "CallAttemptFailure"

  it should "parse metadata" in {
    val res = CallAttemptFailure.buildFailures(failingInSeveralWaysMetadataJson).unsafeRunSync()
    res.size should be(10)
    val (firsts, lasts) = res.partition(_.callFullyQualifiedName == "FailingInSeveralWays.ScatterJobThatWillFailSometimes")
    firsts.size should be(9)
    val last = lasts.head
    firsts.map(_.callFullyQualifiedName).distinct should contain theSameElementsAs
      Vector("FailingInSeveralWays.ScatterJobThatWillFailSometimes")
    firsts.map(_.jobIndex).sorted should contain theSameElementsInOrderAs (41 to 49)
    last.callFullyQualifiedName should be("FailingInSeveralWays.ReadItBackToMe")
    last.workflowId should be("9e83cb0d-0347-431e-8ee6-48f535060845")
  }

  it should "parse empty json as an empty vector" in {
    val res = CallAttemptFailure.buildFailures("{}").unsafeRunSync()
    res should be(empty)
  }

  it should "skipIncompleteAttempts" in {
    val json =
      """{"id": "wf", "calls": {"a.b": [
        |  {"shardIndex":0,"attempt":1,"failures":[{"message":"m"}]},
        |  {"attempt":1}
        |]}}""".stripMargin
    val res = CallAttemptFailure.buildFailures(json).unsafeRunSync()
    res.map(_.message) should be(Vector("m"))
    CallAttemptFailure.buildFailures("\"not json\"").unsafeRunSync() should be(Vector.empty)
  }

  it should "not parse malformed json" in {
    a[ParsingException] should be thrownBy CallAttemptFailure.buildFailures("{").unsafeRunSync()
  }

}

object CallAttemptFailureSpec {
  private val failingInSeveralWaysMetadataJson = {
    val resource = classOf[CallAttemptFailureSpec].getResource("failingInSeveralWaysMetadata.json")
    val path = resource.getPath
    better.files.File(path).contentAsString
  }
}
