package cromwell.util.google

import java.io._
import java.nio.file.Path

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.http.{HttpTransport, InputStreamContent}
import com.google.api.client.json.JsonFactory
import com.google.api.services.storage.Storage
import com.google.api.services.storage.model.Bucket

object GoogleCloudStorage {
  def apply(appName: String, credential: Credential, jsonFactory: JsonFactory, httpTransport: HttpTransport): GoogleCloudStorage = {
    GoogleCloudStorage(new Storage.Builder(httpTransport, jsonFactory, credential).setApplicationName(appName).build())
  }
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

  def uploadObject(gcsPath: GoogleCloudStoragePath, inputStream: InputStream, byteCount: Long) = {
    val mediaContent: InputStreamContent = new InputStreamContent("application/octet-stream", inputStream)
    mediaContent.setLength(byteCount)

    val insertObject = client.objects.insert(gcsPath.bucket, null, mediaContent)
    insertObject.setName(gcsPath.objectName)

    // For small files, let's call setDirectUploadEnabled(true), to
    // reduce the number of HTTP requests made to the server.
    // Here, define small as 2MB:
    if(mediaContent.getLength > 0 && mediaContent.getLength <= 2 * 1000 * 1000) {
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

  def slurpFile(file: GoogleCloudStoragePath, clientSecretsFile: Path): String = {
    new String(downloadObject(file), "UTF-8")
  }
}