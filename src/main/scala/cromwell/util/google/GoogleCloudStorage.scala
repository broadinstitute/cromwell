package cromwell.util.google

import java.io._
import java.math.BigInteger

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.http.{HttpTransport, InputStreamContent}
import com.google.api.client.json.JsonFactory
import com.google.api.client.util.DateTime
import com.google.api.services.storage.Storage
import com.google.api.services.storage.model.Bucket.Owner
import com.google.api.services.storage.model.{Bucket, StorageObject}
import cromwell.util.google.GoogleCloudStorage.GcsBucketInfo

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

  // See comment in uploadObject re small files. Here, define small as 2MB or lower:
  private val smallFileSizeLimit: Long = 2000000

  def uploadObject(gcsPath: GoogleCloudStoragePath, fileContent: String): Unit = {
    val fileBytes = fileContent.getBytes
    val bais = new ByteArrayInputStream(fileBytes)
    uploadObject(gcsPath, bais, fileBytes.length)
  }

  def uploadObject(gcsPath: GoogleCloudStoragePath, inputStream: InputStream, byteCount: Long): Unit = {
    val mediaContent: InputStreamContent = new InputStreamContent("application/octet-stream", inputStream)
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

  def objectSize(gcsPath: GoogleCloudStoragePath): BigInteger = {
    val getObject = client.objects.get(gcsPath.bucket, gcsPath.objectName)
    val storageObject: StorageObject = getObject.execute()
    storageObject.getSize
  }
}