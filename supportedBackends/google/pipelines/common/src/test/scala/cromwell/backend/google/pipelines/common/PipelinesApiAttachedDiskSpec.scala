package cromwell.backend.google.pipelines.common

import common.assertion.CromwellTimeoutSpec
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiAttachedDisk, PipelinesApiEmptyMountedDisk, PipelinesApiWorkingDisk}
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.TryValues
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table

import scala.util.Failure

class PipelinesApiAttachedDiskSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with TryValues {
  val validTable = Table(
    ("unparsed", "parsed"),
    ("/mnt 3 SSD", PipelinesApiEmptyMountedDisk(DiskType.SSD, 3, DefaultPathBuilder.get("/mnt"))),
    ("/mnt/my_path 10 HDD", PipelinesApiEmptyMountedDisk(DiskType.HDD, 10, DefaultPathBuilder.get("/mnt/my_path"))),
    ("local-disk 100 SSD", PipelinesApiWorkingDisk(DiskType.SSD, 100)),
    ("local-disk 100 LOCAL", PipelinesApiWorkingDisk(DiskType.LOCAL, 100))
  )

  it should "parse" in {
    forAll(validTable) { (unparsed, parsed) =>
      PipelinesApiAttachedDisk.parse(unparsed).get shouldEqual parsed
    }
  }

  it should "stringify" in {
    forAll(validTable) { (unparsed, parsed) =>
      parsed.toString shouldEqual unparsed
    }
  }

  val invalidTable = Table(
    "unparsed",
    "local-disk BAD HDD",
    "local-disk 10 BAD",
    "BAD 100 SSD",
    "foobar"
  )

  it should "reject malformed disk mounts" in {
    forAll(invalidTable) { (unparsed) =>
      PipelinesApiAttachedDisk.parse(unparsed) should be(a[Failure[_]])
    }
  }
}
