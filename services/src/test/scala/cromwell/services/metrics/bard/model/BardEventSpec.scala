package cromwell.services.metrics.bard.model
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

class BardEventSpec extends AnyFlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "BardEvent"

  it should "return a map of properties" in {
    val bardEvent = TaskSummaryEvent("done", 0)
    val properties = bardEvent.getProperties
    properties.get("pushToMixpanel") shouldBe false
    properties.get("event") shouldBe "task:summary"
    properties.get("terminationStatus") shouldBe "done"
    properties.get("returnCode") shouldBe 0

  }
}
