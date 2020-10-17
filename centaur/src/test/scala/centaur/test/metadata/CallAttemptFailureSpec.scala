package centaur.test.metadata

import centaur.test.metadata.CallAttemptFailureSpec._
import common.assertion.CromwellTimeoutSpec
import io.circe.ParsingFailure
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class CallAttemptFailureSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "CallAttemptFailure"

  it should "parse metadata" in {
    val res = CallAttemptFailure.buildFailures(failingInSeveralWaysMetadataJson).unsafeRunSync()
    res.size should be(10)
    val firsts = res.take(9)
    val last = res.drop(9).head
    firsts.map(_.callFullyQualifiedName).distinct should contain theSameElementsAs
      Vector("FailingInSeveralWays.ScatterJobThatWillFailSometimes")
    firsts.map(_.jobIndex) should contain theSameElementsInOrderAs (41 to 49)
    last.callFullyQualifiedName should be("FailingInSeveralWays.ReadItBackToMe")
    last.workflowId should be("9e83cb0d-0347-431e-8ee6-48f535060845")
  }

  it should "parse empty json as an empty vector" in {
    val res = CallAttemptFailure.buildFailures("{}").unsafeRunSync()
    res should be(empty)
  }

  it should "not parse malformed json" in {
    val theException: Exception = the[ParsingFailure] thrownBy CallAttemptFailure.buildFailures("{").unsafeRunSync()
    theException should have message "exhausted input"
  }

}

object CallAttemptFailureSpec {
  private val failingInSeveralWaysMetadataJson = {
    val resource = classOf[CallAttemptFailureSpec].getResource("failingInSeveralWaysMetadata.json")
    val path = resource.getPath
    better.files.File(path).contentAsString
  }
}
