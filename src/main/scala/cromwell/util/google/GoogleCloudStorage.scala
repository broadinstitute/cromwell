package cromwell.util.google

import java.io._
import java.math.BigInteger

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.api.client.http.{HttpTransport, InputStreamContent}
import com.google.api.client.json.JsonFactory
import com.google.api.client.util.DateTime
import com.google.api.services.storage.Storage
import com.google.api.services.storage.model.Bucket.Owner
import com.google.api.services.storage.model.{Bucket, StorageObject}
import cromwell.util.google.GoogleCloudStorage.GcsBucketInfo

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

object GoogleCloudStorage {
  def apply(appName: String, credential: Credential, jsonFactory: JsonFactory, httpTransport: HttpTransport): GoogleCloudStorage = {
    GoogleCloudStorage(new Storage.Builder(httpTransport, jsonFactory, credential).setApplicationName(appName).build())
  }

  case class GcsBucketInfo(bucketName: String, location: String, timeCreated: DateTime, owner: Owner)
}

/**
 * Provides a connector which can provide access to Google Cloud Storage
 */
case class GoogleCloudStorage(client: Storage) {
  def listBucket(bucketName: String): GcsBucketInfo = {
    val getBucket = client.buckets().get(bucketName)
    getBucket.setProjection("full")
    val bucket: Bucket = getBucket.execute()

    new GcsBucketInfo(bucketName, bucket.getLocation, bucket.getTimeCreated, bucket.getOwner)
  }

  def listContents(gcsPath: String): Iterable[String] = listContents(GoogleCloudStoragePath(gcsPath))

  def listContents(gcsPath: GoogleCloudStoragePath): Iterable[String] = {
    val listRequest = client.objects().list(gcsPath.bucket)
    listRequest.setPrefix(gcsPath.objectName)

    for {
      listedFile <- listRequest.execute().getItems.asScala
    } yield s"gs://${listedFile.getBucket}/${listedFile.getName}"
  }

  // See comment in uploadObject re small files. Here, define small as 2MB or lower:
  private val smallFileSizeLimit: Long = 2000000

  private def uploadFile(gcsPath: GoogleCloudStoragePath, fileContent: String, contentType: String) = {
    val fileBytes = fileContent.getBytes
    val bais = new ByteArrayInputStream(fileBytes)
    uploadObject(gcsPath, bais, fileBytes.length, contentType)
  }

  def uploadObject(gcsPath: GoogleCloudStoragePath, fileContent: String): Unit = {
    uploadFile(gcsPath, fileContent, "application/octet-stream")
  }

  def uploadJson(gcsPath: GoogleCloudStoragePath, fileContent: String): Unit = {
    uploadFile(gcsPath, fileContent, "application/json")
  }

  def uploadObject(gcsPath: GoogleCloudStoragePath, inputStream: InputStream, byteCount: Long, contentType: String): Unit = {
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

  def deleteObject(gcsPath: GoogleCloudStoragePath): Unit = {
    client.objects.delete(gcsPath.bucket, gcsPath.objectName).execute()
  }

  def downloadObject(gcsPath: GoogleCloudStoragePath): Array[Byte] = {
    val outputStream: ByteArrayOutputStream = new ByteArrayOutputStream()
    val getObject = client.objects.get(gcsPath.bucket, gcsPath.objectName)
    getObject.getMediaHttpDownloader.setDirectDownloadEnabled(true)
    getObject.executeMediaAndDownloadTo(outputStream)

    outputStream.toByteArray
  }

  def slurpFile(file: GoogleCloudStoragePath): String = {
    new String(downloadObject(file), "UTF-8")
  }

  def exists(gcsPath: GoogleCloudStoragePath): Boolean = {
    val getObject = client.objects.get(gcsPath.bucket, gcsPath.objectName)
    Try(getObject.execute) match {
      case Success(_) => true
      case Failure(ex: GoogleJsonResponseException) if ex.getStatusCode == 404 => false
      case Failure(ex) => throw ex
    }
  }

  def objectSize(gcsPath: GoogleCloudStoragePath): BigInteger = {
    val getObject = client.objects.get(gcsPath.bucket, gcsPath.objectName)
    val storageObject: StorageObject = getObject.execute()
    storageObject.getSize
  }
}