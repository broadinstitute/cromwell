package cromwell.languages.util

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class LooksParseableSpec extends AnyFlatSpec with Matchers {
  behavior of "simpleLooksParseable"

  it should "work with `version 1.0` followed by whitespace" in {

    val source =
      """
        |version 1.0              
        |# Programmer warning: make sure your editor does not trim the whitespace in the line above.
        |""".stripMargin

    LanguageFactoryUtil.simpleLooksParseable(List("version 1.0"), List("#"))(source) shouldBe true
  }

  it should "work with `version 1.0` preceded by whitespace" in {

    val source =
      """
        |              version 1.0
        |# Programmer warning: make sure your editor does not trim the whitespace in the line above.
        |""".stripMargin

    LanguageFactoryUtil.simpleLooksParseable(List("version 1.0"), List("#"))(source) shouldBe true
  }

  it should "work with `version 1.0` surrounded by whitespace" in {

    val source =
      """
        |              version 1.0              
        |# Programmer warning: make sure your editor does not trim the whitespace in the line above.
        |""".stripMargin

    LanguageFactoryUtil.simpleLooksParseable(List("version 1.0"), List("#"))(source) shouldBe true
  }

  it should "work with `version 1.0` surrounded by whitespace on all sides" in {

    val source =
      """
        |              
        |              version 1.0              
        |              
        |# Programmer warning: make sure your editor does not trim the whitespace in the line above.
        |""".stripMargin

    LanguageFactoryUtil.simpleLooksParseable(List("version 1.0"), List("#"))(source) shouldBe true
  }

  // This test technically contradicts the spec [0] but there is better hope of receiving a useful
  // error if we try & fail to parse as 1.0, than if we fall back to the server default, `draft-2`
  // > From draft-3 forward, the first line of all WDL files must be a version statement
  // [0] https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#versioning
  it should "work with `version 1.0` surrounded by comments" in {

    val source =
      """
        |# My WDL does a cool thing
        |version 1.0
        |# Here we go...
        |""".stripMargin

    LanguageFactoryUtil.simpleLooksParseable(List("version 1.0"), List("#"))(source) shouldBe true
  }

  it should "reject Chris's idea of a version declaration" in {

    val source =
      """
        |# Thank goodness this WDL is not version 1.1!
        |""".stripMargin

    LanguageFactoryUtil.simpleLooksParseable(List("version 1.0"), List("#"))(source) shouldBe false
  }

}
