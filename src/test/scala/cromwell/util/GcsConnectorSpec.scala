package cromwell.util

import java.nio.file.{Paths, Path}

import org.scalatest._

/**
 * Specification for the GCS Connector
 */
class GcsConnectorSpec extends FlatSpec with Matchers{

  val CLIENT_SECRETS_FILE: Path = Paths.get("/Users/chrisl/client_secrets.json")
  final val TEST_BUCKET = "chrisl-dsde-dev"

  final val TEST_FILE = "testfile"
  final val TEST_FILE_CONTENTS = "test"

  "GcsConnector" should "be able to query GCS buckets for bucket name" in {
    val connector = new GcsConnector("cromwell", CLIENT_SECRETS_FILE)

    val gcsBucketInfo: GcsBucketInfo = connector.listBucket("chrisl-dsde-dev")

    assert(gcsBucketInfo.bucketName == TEST_BUCKET)
  }

  "GcsConnector" should "be able to query GCS buckets for location" in {
    val connector = new GcsConnector("cromwell", CLIENT_SECRETS_FILE)

    val gcsBucketInfo: GcsBucketInfo = connector.listBucket(TEST_BUCKET)

    assert(gcsBucketInfo.location == "US")
  }

  "GcsConnector" should "be able to read file contents" in {
    val connector = new GcsConnector("cromwell", CLIENT_SECRETS_FILE)
    val gcsPath: GoogleCloudStoragePath = new GoogleCloudStoragePath(TEST_BUCKET, TEST_FILE)
    val bytes = connector.downloadObject(gcsPath)

    bytes shouldEqual TEST_FILE_CONTENTS.getBytes
  }

}
