package cromwell.util.google

import org.scalatest.{FlatSpec, Matchers}

/**
 * Spec for GcsPath, in particular 'gs://.../...' URI parsing in constructors.
 */
class GcsPathSpec extends FlatSpec with Matchers {
  final val BUCKET_NAME = "MyBucket"
  final val SIMPLE_OBJECT_NAME = "object.file"
  final val NESTED_OBJECT_NAME = "dir1/dir2/object.file"

  "GcsPath" should "accept Bucket/ObjectName arguments" in {
    val gcsPath: GcsPath = GcsPath(BUCKET_NAME, SIMPLE_OBJECT_NAME)

    gcsPath.bucket shouldEqual BUCKET_NAME
    gcsPath.objectName shouldEqual SIMPLE_OBJECT_NAME
  }

  it should "parse simple gs://Bucket/ObjectName URI strings" in {
    val gcsPath: GcsPath = GcsPath("gs://" + BUCKET_NAME + "/" + SIMPLE_OBJECT_NAME)

    gcsPath.bucket shouldEqual BUCKET_NAME
    gcsPath.objectName shouldEqual SIMPLE_OBJECT_NAME
  }

  it should "parse nested gs://Bucket/dir1/.../dirn/ObjectName URI strings correctly" in {
    val gcsPath: GcsPath = GcsPath("gs://" + BUCKET_NAME + "/" + NESTED_OBJECT_NAME)

    gcsPath.bucket shouldEqual BUCKET_NAME
    gcsPath.objectName shouldEqual NESTED_OBJECT_NAME
  }

  it should "return a Failure if and only if the path cannot be parsed" in {
    val gcsPathTryF = GcsPath.parse("invalid")
    assert(gcsPathTryF.isFailure)

    val gcsPathTryS = GcsPath.parse("gs://valid/path")
    assert(gcsPathTryS.isSuccess)
  }
}
