package cromwell.backend.impl.jes

import java.nio.file.Paths

import com.google.api.services.genomics.model.Disk
import cromwell.backend.impl.jes.io.{DiskType, JesAttachedDisk, JesEmptyMountedDisk, JesWorkingDisk}
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers, TryValues}

import scala.util.Failure

class JesAttachedDiskSpec extends FlatSpec with Matchers with TryValues {
  val validTable = Table(
    ("unparsed", "parsed"),
    ("/mnt 3 SSD", JesEmptyMountedDisk(DiskType.SSD, 3, Paths.get("/mnt"))),
    ("/mnt/my_path 10 HDD", JesEmptyMountedDisk(DiskType.HDD, 10, Paths.get("/mnt/my_path"))),
    ("local-disk 100 SSD", JesWorkingDisk(DiskType.SSD, 100)),
    ("local-disk 100 LOCAL", JesWorkingDisk(DiskType.LOCAL, 100))
  )

  it should "parse" in {
    forAll(validTable) { (unparsed, parsed) =>
      JesAttachedDisk.parse(unparsed).get shouldEqual parsed
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
      JesAttachedDisk.parse(unparsed) should be(a[Failure[_]])
    }
  }

  it should "convert to Google Disk" in {
    val disk = new Disk().setName("d-39de0dbcfb68c8735bd088c62fa061a4")
              .setType(DiskType.SSD.googleTypeName).setAutoDelete(true).setSizeGb(100).setMountPoint("/mnt")
    JesAttachedDisk.parse("/mnt 100 SSD").get.toGoogleDisk shouldEqual disk
  }
}
