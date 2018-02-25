package cromwell.backend.impl.bcs

import org.scalatest.prop.Tables.Table
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.TryValues._

class BcsVpcConfigurationSpec extends BcsTestUtilSpec {
  behavior of s"BcsVpcConfiguration"

  val validVpcTable = Table(
    ("unparsed", "parsed"),
    ("192.168.0.0/16", BcsVpcConfiguration(Some("192.168.0.0/16"), None)),
    ("vpc-xxxx", BcsVpcConfiguration(None, Some("vpc-xxxx"))),
    ("192.168.0.0/16 vpc-xxxx", BcsVpcConfiguration(Some("192.168.0.0/16"), Some("vpc-xxxx"))),
    ("vpc-xxxx 192.168.0.0/16 ", BcsVpcConfiguration(Some("192.168.0.0/16"), Some("vpc-xxxx"))),
  )

  it should "parse correct vpc configuration" in {
    forAll(validVpcTable) { (unparsed, parsed) =>
      BcsVpcConfiguration.parse(unparsed).success.value shouldEqual parsed
    }
  }
}
