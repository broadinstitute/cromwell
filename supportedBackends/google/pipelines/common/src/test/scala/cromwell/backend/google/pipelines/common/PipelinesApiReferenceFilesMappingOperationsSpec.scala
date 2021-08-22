package cromwell.backend.google.pipelines.common

import cats.effect.IO
import com.google.cloud.storage.Storage
import common.assertion.CromwellTimeoutSpec
import cromwell.backend.google.pipelines.common.io.PipelinesApiReferenceFilesDisk
import cromwell.cloudsupport.gcp.auth.MockAuthMode
import cromwell.filesystems.gcs.GcsPathBuilder.ValidFullGcsPath
import org.scalatest.matchers.should.Matchers
import org.scalatest.flatspec.AnyFlatSpecLike

class PipelinesApiReferenceFilesMappingOperationsSpec extends AnyFlatSpecLike with CromwellTimeoutSpec with Matchers {

  private val refFile1Disk1 = "bucketname/dir1/dir2/filename1"
  private val refFile2Disk1 = "bucketname/dir1/dir2/dir3/filename2"
  private val refFile3Disk2 = "bucketname2/dir1/filename"
  private val refFile4Disk2MismatchingChecksum = "bucketname2/dir1/dir2/filename"

  private val disk1 = PipelinesApiReferenceFilesDisk("image_identifier_1", 500)
  private val disk2 = PipelinesApiReferenceFilesDisk("image_identifier_2", 100)

  private val papiReferenceFilesMappingOperationsMockObject: PipelinesApiReferenceFilesMappingOperations =
    new PipelinesApiReferenceFilesMappingOperations {

      override def bulkValidateCrc32cs(gcsClient: Storage, filesWithValidPaths: Map[ReferenceFile, ValidFullGcsPath]): IO[Map[ReferenceFile, Boolean]] =
        IO.pure(filesWithValidPaths.keySet.map(file => (file, file.path != refFile4Disk2MismatchingChecksum)).toMap)
    }

  private val refFileMappingsMock = papiReferenceFilesMappingOperationsMockObject.generateReferenceFilesMapping(
    MockAuthMode("default"),
    List(
      ManifestFile(imageIdentifier = disk1.image, diskSizeGb = disk1.sizeGb, files = List(
        ReferenceFile(path = refFile1Disk1, crc32c = 5),
        ReferenceFile(path = refFile2Disk1, crc32c = 6)
      )),
      ManifestFile(imageIdentifier = disk2.image, diskSizeGb = disk2.sizeGb, files = List(
        ReferenceFile(path = refFile3Disk2, crc32c = 7)
      ))
    )
  )

  it should "correctly figure out which disks have to be mounted based on provided input file paths" in {
    val nonReferenceInputFilePaths = Set("gs://not/a/reference/file")
    val forNonReferenceFile = papiReferenceFilesMappingOperationsMockObject.getReferenceDisksToMount(refFileMappingsMock, nonReferenceInputFilePaths)
    forNonReferenceFile.isEmpty shouldBe true

    val referenceInputFilePathsFrom2Disks = Set(s"gs://$refFile1Disk1", s"gs://$refFile3Disk2")
    val forReferencesFrom2Disks = papiReferenceFilesMappingOperationsMockObject.getReferenceDisksToMount(refFileMappingsMock, referenceInputFilePathsFrom2Disks)
    forReferencesFrom2Disks should contain theSameElementsAs List(disk1, disk2)

    val referenceInputFilePathsFromSingleDisk = Set(s"gs://$refFile1Disk1", s"gs://$refFile2Disk1")
    val forReferencesFromSingleDisk = papiReferenceFilesMappingOperationsMockObject.getReferenceDisksToMount(refFileMappingsMock, referenceInputFilePathsFromSingleDisk)
    forReferencesFromSingleDisk.size shouldBe 1
    forReferencesFromSingleDisk.head shouldBe disk1
  }

  it should "not consider valid a reference file with mismatching checksum" in {
    val mismatchingChecksumReferenceFile = Set(refFile4Disk2MismatchingChecksum)
    val forMismatchingChecksumReferenceFile = papiReferenceFilesMappingOperationsMockObject.getReferenceDisksToMount(refFileMappingsMock, mismatchingChecksumReferenceFile)
    forMismatchingChecksumReferenceFile.isEmpty shouldBe true
  }
}
