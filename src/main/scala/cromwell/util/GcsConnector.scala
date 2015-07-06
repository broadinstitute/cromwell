package cromwell.util

import java.io._

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

/**
 * Provides a connector which can authenticate with Google and provide access to Google Cloud Storage
 *
 * @param applicationName The application name to run under
 * @param secretFile The path to the secret file - must be local and must be JSON
 */
class GcsConnector(applicationName: String, secretFile: String) {

  // Create the data store directory:
  private final val DATA_STORE_DIR: File = new File(System.getProperty("user.home"), ".store/storage_sample")

  // Create the JSON factory:
  private final val JSON_FACTORY: JsonFactory = JacksonFactory.getDefaultInstance

  // Initialise the transport:
  private final val httpTransport: HttpTransport = GoogleNetHttpTransport.newTrustedTransport()

  // Initialise the data store factory:
  private final val dataStoreFactory: FileDataStoreFactory = new FileDataStoreFactory(DATA_STORE_DIR)

  // Authorise:
  private final val credential: Credential = authorise()

  private final val client: Storage = new Storage.Builder(httpTransport, JSON_FACTORY, credential)
    .setApplicationName(this.applicationName)
    .build()

  private def authorise(): Credential =
  {
    val clientSecretsStream = new FileInputStream(secretFile)

    val clientSecrets: GoogleClientSecrets = GoogleClientSecrets.load(
      JSON_FACTORY,
      new InputStreamReader(clientSecretsStream))


    if(clientSecrets.getDetails.getClientId == null || clientSecrets.getDetails.getClientSecret == null)
      {
        throw new Exception("client_secrets not well formed")
      }

    val scopes: java.util.Set[String] = new java.util.HashSet[String]()
    scopes.add(StorageScopes.DEVSTORAGE_READ_WRITE)

    val flow: GoogleAuthorizationCodeFlow = new GoogleAuthorizationCodeFlow.Builder(
      httpTransport,
      JSON_FACTORY,
      clientSecrets,
      scopes
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
    if(mediaContent.getLength > 0 && mediaContent.getLength <= 2 * 1000 * 1000)
      {
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

object GcsConnector {

}