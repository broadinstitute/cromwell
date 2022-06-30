package cromwell.backend.impl.bcs

import org.scalatest.prop.Tables.Table
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.TryValues._

import scala.util.Failure

class BcsDiskSpec extends BcsTestUtilSpec {
  behavior of s"BcsDisk"

  val validDiskTable = Table(
    ("unparsed", "parsed"),
    ("cloud 40", BcsSystemDisk("cloud", 40)),
    ("cloud 200 /home/inputs/", BcsDataDisk("cloud", 200, "/home/inputs/"))
  )

  it should "parse correct disk" in {
    forAll(validDiskTable) { (unparsed, parsed)=>
      BcsDisk.parse(unparsed).success.value shouldEqual(parsed)
    }
  }

  val invalidDiskTable = List(
    "",
    "cloud",
    "40",
    "cloud 40GB",
    "cloud /home/inputs/",
    "cloud /home/inputs 40"
  )

  invalidDiskTable foreach { unparsed =>
    BcsDisk.parse(unparsed) shouldBe a [Failure[_]]
  }
}
