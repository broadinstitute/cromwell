package cromwell.backend.impl.bcs

import org.scalatest.prop.Tables.Table
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.TryValues._


class BcsDockerSpec extends BcsTestUtilSpec {
  behavior of s"BcsDocker"

  val validDockerTable = Table(
    ("unparsed", "parsed"),
    ("ubuntu/latest oss://bcs-reg/ubuntu/", BcsDockerWithPath("ubuntu/latest", "oss://bcs-reg/ubuntu/")),
    ("ubuntu/latest", BcsDockerWithoutPath("ubuntu/latest"))
  )

  it should "parse correct docker configuration" in {
    forAll(validDockerTable) { (unparsed, parsed) =>
      BcsDocker.parse(unparsed).success.value shouldEqual(parsed)
    }
  }
}
