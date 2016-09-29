package lenthall.config

import cats.data.NonEmptyList
import org.scalatest.{FlatSpec, Matchers}


class ConfigValidationExceptionSpec extends FlatSpec with Matchers {

  behavior of "ConfigValidationException"

  it should "Aggregate error messages" in {
    val validationException = ConfigValidationException("You", NonEmptyList.of("Error1",  "Error2"))
    validationException.getMessage shouldBe
      """Invalid You configuration
        |Error1
        |Error2""".stripMargin
  }

  it should "Aggregate error message strings" in {
    val validationException = new ConfigValidationException("You", "Error1 and Error2")
    validationException.getMessage shouldBe
      """Invalid You configuration
        |Error1 and Error2""".stripMargin
  }
}
