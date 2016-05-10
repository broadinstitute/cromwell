package cromwell.filesystems.gcs

import java.io._
import java.nio.file.Paths

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp
import com.google.api.client.googleapis.auth.oauth2.{GoogleAuthorizationCodeFlow, GoogleClientSecrets, GoogleCredential}
import com.google.api.client.googleapis.extensions.java6.auth.oauth2.GooglePromptReceiver
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.client.util.store.FileDataStoreFactory
import com.google.api.services.storage.{Storage, StorageScopes}
import cromwell.filesystems.gcs.GoogleAuthMode.GoogleAuthOptions
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}
import GoogleAuthMode.GcsScopes
import com.typesafe.config.Config

object GoogleAuthMode {

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

  trait GoogleAuthOptions {
    def get(key: String): Try[String]
  }

  val GcsScopes = List(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE
  )
}


sealed trait GoogleAuthMode {

  def credential(options: GoogleAuthOptions): Credential = {
    validateCredentials(buildCredentials(options))
  }

  def assertWorkflowOptions(options: GoogleAuthOptions): Unit = ()

  def name: String

  def requiresAuthFile: Boolean = false

  protected lazy val jsonFactory = JacksonFactory.getDefaultInstance
  protected lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport

  protected lazy val log = LoggerFactory.getLogger(getClass.getSimpleName)

  protected def validateCredentials(credential: Credential) = {
    Try(credential.refreshToken()) match {
      case Failure(ex) => throw new Throwable(s"Google credentials are invalid: ${ex.getMessage}")
      case Success(_) => credential
    }
  }

  protected def buildCredentials(options: GoogleAuthOptions): Credential

  def buildStorage(options: GoogleAuthOptions, config: Config): Storage = {
    new Storage.Builder(
      httpTransport,
      jsonFactory,
      credential(options)).setApplicationName(GoogleConfiguration(config).applicationName).build()
  }
}

final case class ServiceAccountMode(override val name: String, accountId: String, pemPath: String, scopes: List[String] = GcsScopes) extends GoogleAuthMode {
  override protected def buildCredentials(options: GoogleAuthOptions): Credential = {
    val pemFile = new File(pemPath)
    if (!pemFile.exists()) {
      throw new FileNotFoundException(s"PEM file ${Paths.get(pemPath).toAbsolutePath} does not exist")
    }
    new GoogleCredential.Builder().setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setServiceAccountId(accountId)
      .setServiceAccountScopes(scopes.asJava)
      .setServiceAccountPrivateKeyFromPemFile(pemFile)
      .build()
  }
}

final case class UserMode(override val name: String, user: String, secretsFile: String, datastoreDir: String, scopes: List[String] = GcsScopes) extends GoogleAuthMode {
  private def filePathToSecrets(secrets: String, jsonFactory: JsonFactory) = {
    val secretsPath = Paths.get(secrets)
    if(!secretsPath.toFile.canRead) {
      log.warn("Secrets file does not exist or is not readable.")
    }
    val secretStream = new InputStreamReader(new FileInputStream(secretsPath.toFile))

    GoogleClientSecrets.load(jsonFactory, secretStream)
  }

  override protected def buildCredentials(options: GoogleAuthOptions): Credential = {
    val clientSecrets = filePathToSecrets(secretsFile, jsonFactory)
    val dataStore = Paths.get(datastoreDir)
    val dataStoreFactory = new FileDataStoreFactory(dataStore.toFile)
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport,
      jsonFactory,
      clientSecrets,
      scopes.asJava).setDataStoreFactory(dataStoreFactory).build
    new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(user)
  }
}

// It would be goofy to have multiple auths that are application_default, but Cromwell won't prevent it.
final case class ApplicationDefaultMode(override val name: String, scopes: List[String] = GcsScopes) extends GoogleAuthMode {
  override protected def buildCredentials(options: GoogleAuthOptions): Credential = {
    try {
      GoogleCredential.getApplicationDefault().createScoped(scopes.asJava)
    } catch {
      case e: IOException =>
        log.warn("Failed to get application default credentials", e)
        throw e
    }
  }
}

final case class RefreshTokenMode(name: String, clientId: String, clientSecret: String) extends GoogleAuthMode with ClientSecrets {
  import GoogleAuthMode._

  override def requiresAuthFile = true

  /**
    * Throws if the refresh token is not specified.
    */
  override def assertWorkflowOptions(options: GoogleAuthOptions) = getToken(options)

  private def getToken(options: GoogleAuthOptions): String = {
    options.get(RefreshTokenOptionKey).getOrElse(throw new IllegalArgumentException(s"Missing parameters in workflow options: $RefreshTokenOptionKey"))
  }

  override protected def buildCredentials(options: GoogleAuthOptions): Credential = {
    new GoogleCredential.Builder().setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setClientSecrets(clientId, clientSecret)
      .build()
      .setRefreshToken(getToken(options))
  }
}

trait ClientSecrets {
  val clientId: String
  val clientSecret: String
}

final case class SimpleClientSecrets(clientId: String, clientSecret: String) extends ClientSecrets
