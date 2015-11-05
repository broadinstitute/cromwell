package cromwell.util.google

import org.scalatest.{FlatSpec, Matchers}

/**
 * Spec for GcsPath, in particular 'gs://.../...' URI parsing in constructors.
 */
class GcsPathSpec extends FlatSpec with Matchers {
  final val BucketName = "MyBucket"
  final val SimpleObjectName = "object.file"
  final val NestedObjectName = "dir1/dir2/object.file"

  "GcsPath" should "accept Bucket/ObjectName arguments" in {
    val gcsPath: GcsPath = GcsPath(BucketName, SimpleObjectName)

    gcsPath.bucket shouldEqual BucketName
    gcsPath.objectName shouldEqual SimpleObjectName
  }

  it should "parse simple gs://Bucket/ObjectName URI strings" in {
    val gcsPath: GcsPath = GcsPath("gs://" + BucketName + "/" + SimpleObjectName)

    gcsPath.bucket shouldEqual BucketName
    gcsPath.objectName shouldEqual SimpleObjectName
  }

  it should "parse nested gs://Bucket/dir1/.../dirn/ObjectName URI strings correctly" in {
    val gcsPath: GcsPath = GcsPath("gs://" + BucketName + "/" + NestedObjectName)

    gcsPath.bucket shouldEqual BucketName
    gcsPath.objectName shouldEqual NestedObjectName
  }

  it should "return a Failure if and only if the path cannot be parsed" in {
    val gcsPathTryF = GcsPath.parse("invalid")
    assert(gcsPathTryF.isFailure)

    val gcsPathTryS = GcsPath.parse("gs://valid/path")
    assert(gcsPathTryS.isSuccess)
  }
}
