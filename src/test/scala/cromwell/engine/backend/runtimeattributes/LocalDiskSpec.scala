package cromwell.engine.backend.runtimeattributes

import org.scalatest.{TryValues, Matchers, FlatSpec}
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table

class LocalDiskSpec extends FlatSpec with Matchers with TryValues {
  val table = Table(
    ("unparsed", "parsed", "canonical"),
    ("Disk1 3 SSD", LocalDisk("Disk1", DiskType.SSD, 3), "Disk1 3 SSD"),
    ("Disk1 10 HDD", LocalDisk("Disk1", DiskType.HDD, 10), "Disk1 10 HDD"),
    ("foobar LOCAL", LocalDisk("foobar", DiskType.LOCAL), "foobar 10 LOCAL")
  )

  it should "parse" in {
    forAll(table) { (unparsed, parsed, canonical) =>
      LocalDisk.parse(unparsed).get shouldEqual parsed
    }
  }

  it should "stringify" in {
    forAll(table) { (unparsed, parsed, canonical) =>
      parsed.toString shouldEqual canonical
    }
  }
}
