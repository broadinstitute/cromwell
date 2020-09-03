package org.lerch.s3fs

import java.net.URI
import java.util.Properties

import software.amazon.awssdk.services.s3.S3Client

class S3PathSpec extends S3FileSystemUnitSpec {

  val BUCKET_NAME: String = "my-bucket"
  lazy val filesystem: S3FileSystem = {
    val s3Client = S3Client.builder().build()
    val s3fsProvider = new S3FileSystemProvider
    s3fsProvider.createFileSystem(new URI("s3.amazonaws.com"), new Properties, s3Client)
  }

  "An S3Path" should "be created" in {
    val path = new S3Path(filesystem, BUCKET_NAME, "some/interesting/path", "a_file.json" )
    assertResult("my-bucket/some/interesting/path/a_file.json")(path.toString)
  }

  it should " be able to convert to a URI and an encoded string" in {
    val path = new S3Path(filesystem, BUCKET_NAME, "some/interesting#path?or&other", "a_file.json" )
    assertResult(new URI("my-bucket/some/interesting%23path%3For%26other/a_file.json"),
                  "has the Path been properly encoded?")(path.toUri)
    assertResult("my-bucket/some/interesting%23path%3For%26other/a_file.json",
                  "has the Path been URI encoded?")(path.toString)
  }
}
