package common.util

import org.scalatest.{FlatSpec, Matchers}
import common.util.StringUtil.EnhancedToStringable
import common.util.StringUtilSpec.Foo

class StringUtilSpec extends FlatSpec with Matchers {

  it should "correctly truncate a case class with a really long list" in {

    def longLine(i: Int) = s"blah$i" * 100
    val longList = Foo("long long list", (0.until(100)).map(i => longLine(i)).toList)

    longList.toPrettyElidedString(1000) should be("""Foo(
                                                    |  "long long list",
                                                    |  List(
                                                    |    "blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0",
                                                    |    "blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah...""".stripMargin)

  }
}

object StringUtilSpec {
  final case class Foo(bar: String, list: List[String])
}
