package cromwell.api

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import spray.json.JsonParser.ParsingException

class CromwellClientSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
  behavior of "CromwellClient"

  val table = Table(
    ("description", "optionsOption", "refreshTokenOption", "expected"),
    ("ignore bad json when refresh token not provided", Option("{"), None, Option("{")),
    ("ignore bad json when refresh token provided but not used", Option("{"), Option("myToken"),
      Option("{")),
    ("not format json when refresh token key not found", Option("{   }"), Option("myToken"), Option("{   }")),
    ("replace token when found", Option("""{"refresh_token" : "replace_me"}"""), Option("myToken"),
      Option("""{"refresh_token":"myToken"}""")),
  )

  forAll(table) { (description, optionsOption, refreshTokenOption, expected) =>
    it should description in {
      val actual = CromwellClient.insertSecrets(optionsOption, refreshTokenOption)
      actual should be(expected)
    }
  }

  it should "throw an exception when inserting a refresh token into bad json using the token" in {
    val optionsOption = Option("""{"refresh_token" : "replace_me"""")
    val refreshTokenOption = Option("myToken")
    val actual = intercept[ParsingException](CromwellClient.insertSecrets(optionsOption, refreshTokenOption))
    actual.summary should be("""Unexpected end-of-input at input index 31 (line 1, position 32), expected '}'""")
    actual.detail should be(
      """|
         |{"refresh_token" : "replace_me"
         |                               ^
         |""".stripMargin)
  }
}
