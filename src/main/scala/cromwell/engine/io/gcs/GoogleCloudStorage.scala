package cromwell.engine.io.gcs

import java.io._
import java.math.BigInteger

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.api.client.http.InputStreamContent
import com.google.api.client.util.DateTime
import com.google.api.services.storage.Storage
import com.google.api.services.storage.model.Bucket.Owner
import com.google.api.services.storage.model.{Bucket, StorageObject}
import cromwell.engine.PathString
import cromwell.engine.io.IoInterface
import cromwell.engine.workflow.WorkflowOptions
import cromwell.util.TryUtil
import cromwell.util.google.GoogleCredentialFactory

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object GoogleCloudStorage {

  val RefreshTokenOptionKey = "refresh_token"

  lazy val cromwellAuthenticated: Try[GoogleCloudStorage] = {
    for {
      cromwellCredentials <- Try(GoogleCredentialFactory.fromCromwellAuthScheme)
      conf <- GoogleConfiguration.gcloudConf
    } yield apply(conf.appName, cromwellCredentials)
  }

  def userAuthenticated(workflowOptions: WorkflowOptions): Try[GoogleCloudStorage] = for {
    conf <- GoogleConfiguration.gcloudConf
    token <- workflowOptions.get(RefreshTokenOptionKey)
    userCredentials <- GoogleCredentialFactory.fromUserAuthScheme(token)
  } yield apply(conf.appName, userCredentials)

  private def apply(appName: String, credential: Credential): GoogleCloudStorage = {
    GoogleCloudStorage(new Storage.Builder(GoogleCredentialFactory.httpTransport, GoogleCredentialFactory.jsonFactory, credential).setApplicationName(appName).build())
  }

  case class GcsBucketInfo(bucketName: String, location: String, timeCreated: DateTime, owner: Owner)
}

/**
 * Provides a connector which can provide access to Google Cloud Storage
 */
case class GoogleCloudStorage private(client: Storage) extends IoInterface {

  import GcsPath._
  import PathString._

  def isValidPath(path: String) = path.isGcsUrl

  def readFile(path: String): String = {
    new String(downloadObject(path), "UTF-8")
  }

  /**
    * Gets a CRC code from a GCS object. This is a CRC32c checksum which can serve as a hash. It is resistant to
    * composite uploads. See https://cloud.google.com/storage/docs/hashes-etags#_CRC32C for more.
    */
  def getCrc32c(googleCloudStoragePath: GcsPath): String = {
    val obj = client.objects().get(googleCloudStoragePath.bucket, googleCloudStoragePath.objectName).execute()
    obj.getCrc32c
  }

  def hash(path: String) = getCrc32c(GcsPath(path))

  def exists(path: String): Boolean = {
    val getObject = client.objects.get(path.bucket, path.objectName)
    val polling = 500 milliseconds

    def tryExists(retries: Int): Boolean = Try(getObject.execute) match {
      case Success(_) => true
      case Failure(ex: GoogleJsonResponseException) if ex.getStatusCode == 404 && retries > 0 =>
        /* Could not use TryUtil here because it requires a WorkflowLogger, which we can't get form here
         * TODO From a more general perspective we may need to add a logging capability to the IOInterface
         */
        Thread.sleep(polling.toMillis)
        tryExists(retries - 1)
      case Failure(ex: GoogleJsonResponseException) if ex.getStatusCode == 404 => false
      case Failure(ex) => throw ex
    }

    tryExists(3)
  }

  def listContents(path: String): Iterable[String] = {
    val listRequest = client.objects().list(path.bucket)
    listRequest.setPrefix(path.objectName)

    for {
      listedFile <- listRequest.execute().getItems.asScala
    } yield s"gs://${listedFile.getBucket}/${listedFile.getName}"
  }

  def writeFile(path: String, content: String): Unit = uploadObject(path, content)

  /**
    * Copy file from one GCS path to another
    *
    * @param from - source GCS path (must exist and point to a file)
    * @param to - destination GCS path
    * @return a Try[StorageObject] which is a result of the call to Storage.Objects.copy()
    */
  def copy(from: GcsPath, to: GcsPath): Try[StorageObject] = Try {
    val storageObject = client.objects.get(from.bucket, from.objectName).execute
    client.objects.copy(from.bucket, from.objectName, to.bucket, to.objectName, storageObject).execute
  }

  /**
    * Copy file from one GCS path to another
    *
    * @param from - source GCS path (must exist and point to a file)
    * @param to - destination GCS path
    * @return a Try[StorageObject] which is a result of the call to Storage.Objects.copy()
    */
  def copy(from: NioGcsPath, to: NioGcsPath): Try[StorageObject] = Try {
    val storageObject = client.objects.get(from.bucket, from.objectName).execute
    client.objects.copy(from.bucket, from.objectName, to.bucket, to.objectName, storageObject).execute
  }

  /**
    * Copy files with prefix `from` to files with prefix `to`
    *
    * For example, calling copyPrefix("gs://bucket/a/", "gs://other/prefix/"):
    *
    * 1) listContents(from) = Seq("gs://bucket/a/b.txt", "gs://bucket/a/b/c.txt")
    * 2) The following copies will take place:
    *     a) gs://bucket/a/b.txt -> gs://other/prefix/a/b.txt
    *     b) gs://bucket/a/b/c.txt -> gs://other/prefix/a/b/c.txt
    *
    * @param from - GCS URL prefix to copy files from
    * @param to - GCS URL prefix to copy files to
    */
  def copy(from: String, to: String): Unit = {
    Try(listContents(from)) flatMap { paths =>
      TryUtil.sequence(paths map {path =>
        val source = GcsPath(path)
        val dest = GcsPath(path.replaceAll(s"^$from", to))
        copy(source, dest)
      } toSeq)
    }
  }

  def move(from: NioGcsPath, to: NioGcsPath): Unit = {
    val polling = 500 milliseconds
    def moveInner = {
      val storageObject = client.objects.get(from.bucket, from.objectName).execute
      client.objects.rewrite(from.bucket, from.objectName, to.bucket, to.objectName, storageObject).execute
    }

    def tryMove(retries: Int): Unit = Try(moveInner) match {
      case Success(_) =>
      case Failure(ex: GoogleJsonResponseException) if ex.getStatusCode == 404 && retries > 0 =>
        /* Could not use TryUtil here because it requires a WorkflowLogger, which we can't get form here
         * TODO From a more general perspective we may need to add a logging capability to the IOInterface
         */
        Thread.sleep(polling.toMillis)
        tryMove(retries - 1)
      case Failure(ex: GoogleJsonResponseException) if ex.getStatusCode == 404 =>
      case Failure(ex) => throw ex
    }

    tryMove(3)
  }

  //TODO: improve to honor pattern ?
  def glob(path: String, pattern: String): Seq[String] = listContents(path).toSeq

  def listBucket(bucketName: String): GcsBucketInfo = {
    val getBucket = client.buckets().get(bucketName)
    getBucket.setProjection("full")
    val bucket: Bucket = getBucket.execute()

    new GcsBucketInfo(bucketName, bucket.getLocation, bucket.getTimeCreated, bucket.getOwner)
  }

  // See comment in uploadObject re small files. Here, define small as 2MB or lower:
  private val smallFileSizeLimit: Long = 2000000

  private def uploadFile(gcsPath: GcsPath, fileContent: String, contentType: String) = {
    val fileBytes = fileContent.getBytes
    val bais = new ByteArrayInputStream(fileBytes)
    uploadObject(gcsPath, bais, fileBytes.length, contentType)
  }

  def uploadObject(gcsPath: GcsPath, fileContent: String): Unit = {
    uploadFile(gcsPath, fileContent, "application/octet-stream")
  }

  def uploadJson(gcsPath: GcsPath, fileContent: String): Unit = {
    uploadFile(gcsPath, fileContent, "application/json")
  }

  def uploadObject(gcsPath: GcsPath, inputStream: InputStream, byteCount: Long, contentType: String): Unit = {
    val mediaContent: InputStreamContent = new InputStreamContent(contentType, inputStream)
    mediaContent.setLength(byteCount)

    val insertObject = client.objects.insert(gcsPath.bucket, null, mediaContent)
    insertObject.setName(gcsPath.objectName)

    // For small files, let's call setDirectUploadEnabled(true), to
    // reduce the number of HTTP requests made to the server.
    if (mediaContent.getLength > 0 && mediaContent.getLength <= smallFileSizeLimit) {
      insertObject.getMediaHttpUploader.setDirectUploadEnabled(true)
    }

    insertObject.execute()
  }

  def deleteObject(gcsPath: GcsPath): Unit = client.objects.delete(gcsPath.bucket, gcsPath.objectName).execute()

  def deleteObject(gcsPath: NioGcsPath): Unit = client.objects.delete(gcsPath.bucket, gcsPath.objectName).execute()

  def downloadObject(gcsPath: GcsPath): Array[Byte] = {
    val outputStream: ByteArrayOutputStream = new ByteArrayOutputStream()
    val getObject = client.objects.get(gcsPath.bucket, gcsPath.objectName)
    getObject.getMediaHttpDownloader.setDirectDownloadEnabled(true)
    getObject.executeMediaAndDownloadTo(outputStream)

    outputStream.toByteArray
  }

  override def size(path: String): Long = objectSize(GcsPath(path)).longValue()

  def objectSize(gcsPath: GcsPath): BigInteger = {
    val getObject = client.objects.get(gcsPath.bucket, gcsPath.objectName)
    val storageObject: StorageObject = getObject.execute()
    storageObject.getSize
  }
}
