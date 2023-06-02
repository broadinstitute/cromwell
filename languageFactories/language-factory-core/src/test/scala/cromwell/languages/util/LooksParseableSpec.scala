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

  it should "reject Chris's idea of a version declaration" in {

    val source =
      """
        |# Thank goodness this WDL is not version 1.1!
        |""".stripMargin

    LanguageFactoryUtil.simpleLooksParseable(List("version 1.0"), List("#"))(source) shouldBe false
  }

}
