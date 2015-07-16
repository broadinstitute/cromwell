package cromwell.util

import cromwell.util.google.GoogleCloudStoragePath
import org.scalatest.{FlatSpec, Matchers}

/**
 * Spec for GoogleCloudStoragePath, in particular 'gs://.../...' URI parsing in constructors.
 */
class GoogleCloudStoragePathSpec extends FlatSpec with Matchers{

  final val BUCKET_NAME = "MyBucket"
  final val SIMPLE_OBJECT_NAME = "object.file"
  final val NESTED_OBJECT_NAME = "dir1/dir2/object.file"

  "GoogleCloudStoragePath" should "accept Bucket/ObjectName arguments" in {
    val gcsPath: GoogleCloudStoragePath = GoogleCloudStoragePath(BUCKET_NAME, SIMPLE_OBJECT_NAME)

    gcsPath.bucket shouldEqual BUCKET_NAME
    gcsPath.objectName shouldEqual SIMPLE_OBJECT_NAME
  }

  "GoogleCloudStoragePath" should "parse simple gs://Bucket/ObjectName URI strings" in {
    val gcsPath: GoogleCloudStoragePath = GoogleCloudStoragePath("gs://" + BUCKET_NAME + "/" + SIMPLE_OBJECT_NAME)

    gcsPath.bucket shouldEqual BUCKET_NAME
    gcsPath.objectName shouldEqual SIMPLE_OBJECT_NAME
  }

  "GoogleCloudStoragePath" should "parse nested gs://Bucket/dir1/.../dirn/ObjectName URI strings correctly" in {
    val gcsPath: GoogleCloudStoragePath = GoogleCloudStoragePath("gs://" + BUCKET_NAME + "/" + NESTED_OBJECT_NAME)

    gcsPath.bucket shouldEqual BUCKET_NAME
    gcsPath.objectName shouldEqual NESTED_OBJECT_NAME
  }

  "GoogleCloudStoragePath's tryParse" should "return a Failure if and only if the path cannot be parsed" in {
    val gcsPathTryF = GoogleCloudStoragePath.parse("invalid")
    assert(gcsPathTryF.isFailure)

    val gcsPathTryS = GoogleCloudStoragePath.parse("gs://valid/path")
    assert(gcsPathTryS.isSuccess)
  }
}
