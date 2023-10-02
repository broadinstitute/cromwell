package cromwell.backend.google.pipelines.v2beta

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class PipelinesUtilityConversionsSpec extends AnyFlatSpec with Matchers {
  behavior of "PipelinesUtilityConversions"

  it should "not modify strings that contain only ascii characters" in {
    val input = "hi there!?"
    PipelinesUtilityConversions.cleanUtf8mb4(input) shouldBe input
  }

  it should "not modify strings with 3-byte unicode characters" in {
    val input = "Here is my non-ascii character: \u1234 Do you like it?"
    PipelinesUtilityConversions.cleanUtf8mb4(input) shouldBe input
  }

  it should "replace 4-byte unicode characters" in {
    val cry = new String(Character.toChars(Integer.parseInt("1F62D", 16)))
    val barf = new String(Character.toChars(Integer.parseInt("1F92E", 16)))
    val input = s"When I try to put an emoji in the database it $barf and then I $cry"
    val cleaned = "When I try to put an emoji in the database it \uFFFD and then I \uFFFD"
    PipelinesUtilityConversions.cleanUtf8mb4(input) shouldBe cleaned
  }
}
