package cromwell.util

import java.io._
import java.nio.file.Path
import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.extensions.java6.auth.oauth2.{AuthorizationCodeInstalledApp, VerificationCodeReceiver}
import com.google.api.client.googleapis.auth.oauth2.{GoogleAuthorizationCodeFlow, GoogleClientSecrets}
import com.google.api.client.googleapis.extensions.java6.auth.oauth2.GooglePromptReceiver
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.http.{InputStreamContent, HttpTransport}
import com.google.api.client.json.JsonFactory
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.client.util.store.FileDataStoreFactory
import com.google.api.services.storage.model.Bucket
import com.google.api.services.storage.{Storage, StorageScopes}
import scala.collection.JavaConverters._

/**
 * Provides a connector which can authenticate with Google and provide access to Google Cloud Storage
 *
 * @param applicationName The application name to run under
 * @param secretFile The path to the secret file - must be local and must be JSON
 */
case class GcsConnector(applicationName: String, secretFile: Path) {

  // Create the data store directory:
  private val dataStoreDirectory: File = new File(System.getProperty("user.home"), ".store/storage_sample")

  // Create the JSON factory:
  private val jsonFactory: JsonFactory = JacksonFactory.getDefaultInstance

  // Initialise the transport:
  private val httpTransport: HttpTransport = GoogleNetHttpTransport.newTrustedTransport()

  // Initialise the data store factory:
  private val dataStoreFactory: FileDataStoreFactory = new FileDataStoreFactory(dataStoreDirectory)

  // Get cred:
  private val credential: Credential = authorize()

  private val client: Storage = new Storage.Builder(httpTransport, jsonFactory, credential)
    .setApplicationName(this.applicationName)
    .build()

  private def authorize(): Credential =
  {
    val clientSecretsStream = new FileInputStream(secretFile.toString)

    val clientSecrets: GoogleClientSecrets = GoogleClientSecrets.load(
      jsonFactory,
      new InputStreamReader(clientSecretsStream))

    if(Option(clientSecrets.getDetails.getClientId).isEmpty || Option(clientSecrets.getDetails.getClientSecret).isEmpty)
      {
        throw new IllegalArgumentException("")
      }

    val scopes: Set[String] = Set(StorageScopes.DEVSTORAGE_READ_WRITE)

    val flow: GoogleAuthorizationCodeFlow = new GoogleAuthorizationCodeFlow.Builder(
      httpTransport,
      jsonFactory,
      clientSecrets,
      scopes.asJava
    )
      .setDataStoreFactory(dataStoreFactory)
      .build()

    val receiver: VerificationCodeReceiver = new GooglePromptReceiver()

    new AuthorizationCodeInstalledApp(flow, receiver).authorize("user")
  }

  def listBucket(bucketName: String): GcsBucketInfo =
  {
    val getBucket = client.buckets().get(bucketName)
    getBucket.setProjection("full")
    val bucket: Bucket = getBucket.execute()

    new GcsBucketInfo(bucketName, bucket.getLocation, bucket.getTimeCreated, bucket.getOwner)
  }

  def uploadObject(gcsPath: GoogleCloudStoragePath, inputStream: InputStream, byteCount: Long) =
  {
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

  def downloadObject(gcsPath: GoogleCloudStoragePath): Array[Byte] =
  {
    val outputStream: ByteArrayOutputStream = new ByteArrayOutputStream()
    val getObject = client.objects.get(gcsPath.bucket, gcsPath.objectName)
    getObject.getMediaHttpDownloader.setDirectDownloadEnabled(true)
    getObject.executeMediaAndDownloadTo(outputStream)

    outputStream.toByteArray
  }
}