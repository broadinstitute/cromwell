package lenthall.config

import org.scalatest.{Matchers, FlatSpec}

import scalaz.NonEmptyList

class ConfigValidationExceptionSpec extends FlatSpec with Matchers {

  "ConfigValidationException" should "Aggregate error messages" in {
    val validationException = ConfigValidationException("You", NonEmptyList("Error1",  "Error2"))
    validationException.getMessage shouldBe
      """Invalid You configuration
        |Error1
        |Error2""".stripMargin
  }

}
