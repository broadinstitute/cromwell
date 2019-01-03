package cromwell.backend.google.pipelines.v1alpha2

import com.google.api.services.genomics.model.Disk
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiAttachedDisk}
import cromwell.backend.google.pipelines.v1alpha2.PipelinesConversions._
import org.scalatest.{FlatSpec, Matchers, TryValues}

class PipelinesAttachedDiskSpec extends FlatSpec with Matchers with TryValues {
  it should "convert to Google Disk" in {
    val disk = new Disk().setName("d-39de0dbcfb68c8735bd088c62fa061a4")
              .setType(DiskType.SSD.googleTypeName).setAutoDelete(true).setSizeGb(100).setMountPoint("/mnt")
    PipelinesApiAttachedDisk.parse("/mnt 100 SSD").get.toGoogleDisk shouldEqual disk
  }
}
