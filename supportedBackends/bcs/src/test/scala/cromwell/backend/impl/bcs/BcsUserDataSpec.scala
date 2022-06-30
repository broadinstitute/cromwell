package cromwell.backend.impl.bcs

import org.scalatest.TryValues._

class BcsUserDataSpec extends BcsTestUtilSpec {
  behavior of s"BcsUserDataSpec"

  it should "work for right user data" in {
    val key = "key"
    val value = "value"

    val userData = BcsUserData.parse(s"$key $value")
    userData.success.value.key shouldEqual key
    userData.success.value.value shouldEqual value

    BcsUserData.parse(s"$key$value").isFailure shouldBe true
  }
}
