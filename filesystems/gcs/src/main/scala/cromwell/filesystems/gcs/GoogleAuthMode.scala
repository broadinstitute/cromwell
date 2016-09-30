package cromwell.filesystems.gcs

import java.io.{FileNotFoundException, IOException, InputStreamReader}
import java.nio.file.{Files, Paths}

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp
import com.google.api.client.googleapis.auth.oauth2.{GoogleAuthorizationCodeFlow, GoogleClientSecrets, GoogleCredential}
import com.google.api.client.googleapis.extensions.java6.auth.oauth2.GooglePromptReceiver
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.client.util.store.FileDataStoreFactory
import com.google.api.services.storage.{Storage, StorageScopes}
import cromwell.filesystems.gcs.GoogleAuthMode.{GcsScopes, GoogleAuthOptions}
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

object GoogleAuthMode {

  lazy val jsonFactory = JacksonFactory.getDefaultInstance
  lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport
  val RefreshTokenOptionKey = "refresh_token"

  /**
    * Before it returns the raw credential, checks if the token will expire within 60 seconds.
    *
    * TODO: Needs more design / testing around thread safety.
    * For example, the credential returned is mutable, and may be modified by another thread.
    *
    * Most Google clients have the ability to refresh tokens automatically, as they use the standard Google
    * HttpTransport that automatically triggers credential refreshing via Credential.handleResponse. Since Cromwell
    * contacts https://gcr.io directly via HTTP requests using spray-client, we need to keep the token fresh ourselves.
    *
    * @see Credential#handleResponse(HttpRequest, HttpResponse, boolean)
    */
  implicit class EnhancedCredentials(val credential: Credential) extends AnyVal {
    def freshCredential: Try[Credential] = {
      val stillValid = Option(credential.getExpiresInSeconds).exists(_ > 60)
      if (stillValid || credential.refreshToken()) {
        Success(credential)
      } else {
        Failure(new Exception("Unable to refresh token"))
      }
    }
  }

  def buildStorage(credential: Credential, applicationName: String) = {
    new Storage.Builder(
      httpTransport,
      jsonFactory,
      credential).setApplicationName(applicationName).build()
  }

  trait GoogleAuthOptions {
    def get(key: String): Try[String]
  }

  val GcsScopes = List(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE
  )
}


sealed trait GoogleAuthMode {
  def credential(options: GoogleAuthOptions): Credential

  def assertWorkflowOptions(options: GoogleAuthOptions): Unit = ()

  def name: String

  def requiresAuthFile: Boolean = false

  protected lazy val log = LoggerFactory.getLogger(getClass.getSimpleName)

  protected def validateCredentials(credential: Credential) = {
    Try(credential.refreshToken()) match {
      case Failure(ex) => throw new RuntimeException(s"Google credentials are invalid: ${ex.getMessage}")
      case Success(_) => credential
    }
  }

  def buildStorage(options: GoogleAuthOptions, applicationName: String): Storage = {
    GoogleAuthMode.buildStorage(credential(options), applicationName)
  }
}

final case class ServiceAccountMode(override val name: String, accountId: String, pemPath: String, scopes: List[String] = GcsScopes) extends GoogleAuthMode {
  import GoogleAuthMode._

  private lazy val credentials: Credential = {
    val pemFile = Paths.get(pemPath).toAbsolutePath
    if (!Files.exists(pemFile)) {
      throw new FileNotFoundException(s"PEM file $pemFile does not exist")
    }
    validateCredentials(
      new GoogleCredential.Builder().setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setServiceAccountId(accountId)
      .setServiceAccountScopes(scopes.asJava)
      .setServiceAccountPrivateKeyFromPemFile(pemFile.toFile)
      .build()
    )
  }

  override def credential(options: GoogleAuthOptions) = credentials
}

final case class UserMode(override val name: String, user: String, secretsFile: String, datastoreDir: String, scopes: List[String] = GcsScopes) extends GoogleAuthMode {
  import GoogleAuthMode._

  private def filePathToSecrets(secrets: String, jsonFactory: JsonFactory) = {
    val secretsPath = Paths.get(secrets).toAbsolutePath
    if(!Files.isReadable(secretsPath)) {
      log.warn("Secrets file does not exist or is not readable.")
    }
    val secretStream = new InputStreamReader(Files.newInputStream(secretsPath))

    GoogleClientSecrets.load(jsonFactory, secretStream)
  }

  private lazy val credentials: Credential = {
    val clientSecrets = filePathToSecrets(secretsFile, jsonFactory)
    val dataStore = Paths.get(datastoreDir).toAbsolutePath
    val dataStoreFactory = new FileDataStoreFactory(dataStore.toFile)
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport,
      jsonFactory,
      clientSecrets,
      scopes.asJava).setDataStoreFactory(dataStoreFactory).build
    validateCredentials(new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(user))
  }

  override def credential(options: GoogleAuthOptions) = credentials
}

// It would be goofy to have multiple auths that are application_default, but Cromwell won't prevent it.
final case class ApplicationDefaultMode(override val name: String, scopes: List[String] = GcsScopes) extends GoogleAuthMode {
  import GoogleAuthMode._

  private lazy val credentials: Credential = {
    try {
      validateCredentials(GoogleCredential.getApplicationDefault().createScoped(scopes.asJava))
    } catch {
      case e: IOException =>
        log.warn("Failed to get application default credentials", e)
        throw e
    }
  }

  override def credential(options: GoogleAuthOptions) = credentials
}

final case class RefreshTokenMode(name: String, clientId: String, clientSecret: String) extends GoogleAuthMode with ClientSecrets {
  import GoogleAuthMode._

  override def requiresAuthFile = true

  /**
    * Throws if the refresh token is not specified.
    */
  override def assertWorkflowOptions(options: GoogleAuthOptions): Unit = { getToken(options); () }

  private def getToken(options: GoogleAuthOptions): String = {
    options.get(RefreshTokenOptionKey).getOrElse(throw new IllegalArgumentException(s"Missing parameters in workflow options: $RefreshTokenOptionKey"))
  }

  override def credential(options: GoogleAuthOptions): Credential = {
    validateCredentials(
      new GoogleCredential.Builder().setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setClientSecrets(clientId, clientSecret)
      .build()
      .setRefreshToken(getToken(options))
    )
  }
}

trait ClientSecrets {
  val clientId: String
  val clientSecret: String
}

final case class SimpleClientSecrets(clientId: String, clientSecret: String) extends ClientSecrets
