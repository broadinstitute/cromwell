package centaur.test

import com.azure.storage.blob.BlobContainerClient
import com.google.cloud.storage.Storage.BlobListOption
import com.google.cloud.storage.{Blob, Storage}
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.model.ListObjectsRequest

import scala.jdk.CollectionConverters._
import scala.language.implicitConversions

trait ObjectCounter[A] {
  def parsePath(regex: String): String => Path = { fullPath =>
    if (fullPath.matches(regex)) {
      val prefixLength = 5
      val bucketAndDashes = fullPath.drop(prefixLength).split("/", 2)
      val bucket = bucketAndDashes.head
      val directory = bucketAndDashes.tail.mkString
      Path(bucket, directory)
    } else throw IllegalPathException()
  }

  def countObjectsAtPath(client: A): Path => Int
}

object ObjectCounterInstances {

  implicit val awsS3ObjectCounter: ObjectCounter[S3Client] = (amazonS3: S3Client) => {
    implicit def boolToInt(b: Boolean) = if (b) 1 else 0

    def isFileExists(s3: S3Client, bucket: String, prefix: String): Boolean =
      !s3.listObjects(ListObjectsRequest.builder().bucket(bucket).prefix(prefix).build()).contents().isEmpty

    AwsS3Path => isFileExists(amazonS3, AwsS3Path.bucket, AwsS3Path.directory).toInt
  }

  implicit val gcsObjectCounter: ObjectCounter[Storage] = (storage: Storage) => {
    val listObjectsAtPath: Path => Iterable[Blob] = g =>
      storage.list(g.bucket, BlobListOption.prefix(g.directory)).iterateAll.asScala
    listObjectsAtPath(_).size
  }

  implicit val blobObjectCounter: ObjectCounter[BlobContainerClient] = (containerClient: BlobContainerClient) => {
    val pathToInt: Path => Int = providedPath => {
      // Our path parsing is somewhat GCP centric. Convert to a blob path starting from the container root.
      def pathToBlobPath(parsedPath: Path): String =
        (Option(parsedPath.bucket), Option(parsedPath.directory)) match {
          case (None, _) => ""
          case (Some(_), None) => parsedPath.bucket
          case (Some(_), Some(_)) => parsedPath.bucket + "/" + parsedPath.directory
        }

      val fullPath = pathToBlobPath(providedPath)
      val blobsInFolder = containerClient.listBlobsByHierarchy(fullPath)
      // if something "isPrefix", it's a directory. Otherwise, its a file. We just want to count files.
      blobsInFolder.asScala.count(!_.isPrefix)
    }
    pathToInt(_)
  }
}

object ObjectCounterSyntax {

  implicit class ObjectCounterSyntax[A](client: A) {
    def countObjects(regex: String)(implicit c: ObjectCounter[A]): String => Int =
      c.parsePath(regex) andThen c.countObjectsAtPath(client)
  }

}

final case class Path(bucket: String, directory: String)

final case class IllegalPathException(private val message: String = "This is not valid path") extends Exception(message)
