package common.util

import common.assertion.CromwellTimeoutSpec
import common.util.StringUtil.EnhancedToStringable
import common.util.StringUtilSpec._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class StringUtilSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  it should "correctly truncate a case class with a really long list" in {
    val fooOfBars = Foo("long long list", 0.until(100).toList.map(i => new Bar(i)))

    // If we try the naive toString, we get an exception when we toString the later elements:
    assertThrows[BarException](fooOfBars.toString)

    // With the elided string, we stop processing early and are able to produce a nice, short string without ever
    // touching the later elements:
    fooOfBars.toPrettyElidedString(1000) should be("""Foo(
                                                    |  "long long list",
                                                    |  List(
                                                    |    "blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0",
                                                    |    "blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah...""".stripMargin)

  }
}

object StringUtilSpec {
  final case class Foo(bar: String, list: List[Bar])

  final class Bar(index: Int) {
    private def longLine(i: Int) = '"' +  s"blah$i" * 100 + '"'
    override def toString: String = if (index < 2) {
      longLine(index)
    } else {
      throw BarException(s"Don't look at index $index!")
    }
  }

  final case class BarException(msg: String) extends Exception {
    override def getMessage: String = msg
  }
}
