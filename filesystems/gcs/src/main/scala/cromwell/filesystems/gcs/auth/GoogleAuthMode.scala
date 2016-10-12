package cromwell.filesystems.gcs.auth

import java.io.{FileNotFoundException, InputStreamReader}
import java.nio.file.Paths

import better.files._
import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp
import com.google.api.client.googleapis.auth.oauth2.{GoogleAuthorizationCodeFlow, GoogleClientSecrets, GoogleCredential}
import com.google.api.client.googleapis.extensions.java6.auth.oauth2.GooglePromptReceiver
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.googleapis.testing.auth.oauth2.MockGoogleCredential
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.client.util.store.FileDataStoreFactory
import com.google.api.services.storage.StorageScopes
import com.google.auth.oauth2.{ClientId, ServiceAccountCredentials}
import com.google.cloud.AuthCredentials
import cromwell.core.WorkflowOptions
import cromwell.filesystems.gcs.auth.GoogleAuthMode._
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

object GoogleAuthMode {

  lazy val jsonFactory = JacksonFactory.getDefaultInstance
  lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport

  val RefreshTokenOptionKey = "refresh_token"
  val GcsScopes = List(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE
  ).asJava

  def checkReadable(file: File) = {
    if (!file.isReadable) throw new FileNotFoundException(s"File $file does not exist or is not readable")
  }

  case object NoAuthMode extends GoogleAuthMode {
    override def name = "no_auth"

    override def authCredentials(options: WorkflowOptions): AuthCredentials = AuthCredentials.noAuth()
    override def credential(options: WorkflowOptions): Credential = new MockGoogleCredential.Builder().build()
  }
}


sealed trait GoogleAuthMode {
  protected lazy val log = LoggerFactory.getLogger(getClass.getSimpleName)

  /**
    * Validate the auth mode against provided options
    */
  def validate(options: WorkflowOptions): Unit = {()}

  def name: String
  // Create an AuthCredentials object from the google-cloud library (https://github.com/GoogleCloudPlatform/google-cloud-java using https://github.com/google/google-auth-library-java under the hood)
  def authCredentials(options: WorkflowOptions): AuthCredentials
  // Create a Credential object from the google.api.client.auth library (https://github.com/google/google-api-java-client)
  def credential(options: WorkflowOptions): Credential

  def requiresAuthFile: Boolean = false

  protected def validateAuthCredentials(authCredentials: AuthCredentials, scopes: java.util.Collection[String]): AuthCredentials = validate(authCredentials, authCredentials.credentials().createScoped(scopes).refresh)

  protected def validateCredential(credential: Credential) = validate(credential, credential.refreshToken)

  private def validate[T](credential: T, validation: () => Any): T = {
    Try(validation()) match {
      case Failure(ex) => throw new RuntimeException(s"Google credentials are invalid: ${ex.getMessage}")
      case Success(_) => credential
    }
  }
}

final case class ServiceAccountMode(override val name: String,
                                    accountId: String,
                                    pemPath: String,
                                    scopes: java.util.List[String]) extends GoogleAuthMode {
  private val pemFile = File(pemPath)
  checkReadable(pemFile)

  private lazy val _authCredentials: AuthCredentials = {
    val saCredentials = ServiceAccountCredentials.fromPkcs8(accountId, accountId, pemFile.contentAsString, null, scopes)
    validateAuthCredentials(AuthCredentials.createFor(saCredentials.getClientId, saCredentials.getPrivateKey), scopes)
  }

  private lazy val _credential: Credential = {
    validateCredential(
      new GoogleCredential.Builder().setTransport(httpTransport)
        .setJsonFactory(jsonFactory)
        .setServiceAccountId(accountId)
        .setServiceAccountScopes(scopes)
        .setServiceAccountPrivateKeyFromPemFile(pemFile.toJava)
        .build()
    )
  }

  override def authCredentials(options: WorkflowOptions) = _authCredentials

  override def credential(options: WorkflowOptions): Credential = _credential
}

final case class UserMode(override val name: String,
                          user: String,
                          val secretsPath: String,
                          datastoreDir: String,
                          scopes: java.util.List[String]) extends GoogleAuthMode {

  private lazy val secrets = {
    val secretsFile = File(secretsPath)
    checkReadable(secretsFile)

    val secretStream = new InputStreamReader(secretsFile.newInputStream)

    GoogleClientSecrets.load(jsonFactory, secretStream)
  }

  private lazy val _credential: Credential = {
    val dataStore = Paths.get(datastoreDir).toAbsolutePath
    val dataStoreFactory = new FileDataStoreFactory(dataStore.toFile)
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport, jsonFactory, secrets, scopes).setDataStoreFactory(dataStoreFactory).build
    validateCredential(new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(user))
  }

  private lazy val _authCredentials: AuthCredentials = {
    new RefreshableOAuth2Credentials(_credential.getRefreshToken, new ClientId(secrets.getDetails.getClientId, secrets.getDetails.getClientSecret))
  }

  override def credential(options: WorkflowOptions) = _credential

  override def authCredentials(options: WorkflowOptions) = _authCredentials
}

private object ApplicationDefault {
  private [auth] lazy val _AuthCredentials = AuthCredentials.createApplicationDefaults()
  private [auth] lazy val _Credential: Credential = GoogleCredential.getApplicationDefault()
}

final case class ApplicationDefaultMode(name: String) extends GoogleAuthMode {
  override def authCredentials(options: WorkflowOptions) = ApplicationDefault._AuthCredentials
  override def credential(options: WorkflowOptions) = ApplicationDefault._Credential
}

final case class RefreshTokenMode(name: String,
                                  clientId: String,
                                  clientSecret: String,
                                  scopes: java.util.List[String]) extends GoogleAuthMode with ClientSecrets {
  import GoogleAuthMode._
  override def requiresAuthFile = true

  private def extractRefreshToken(options: WorkflowOptions): String = {
    options.get(RefreshTokenOptionKey) getOrElse {
      throw new IllegalArgumentException(s"Missing parameters in workflow options: $RefreshTokenOptionKey")
    }
  }

  override def validate(options: WorkflowOptions) = {
    extractRefreshToken(options)

    ()
  }

  override def authCredentials(options: WorkflowOptions): AuthCredentials = {
    val refreshToken = extractRefreshToken(options)
    validateAuthCredentials(new RefreshableOAuth2Credentials(refreshToken, new ClientId(clientId, clientSecret)), scopes)
  }

  override def credential(options: WorkflowOptions): Credential = {
    val refreshToken = extractRefreshToken(options)
    validateCredential(
      new GoogleCredential.Builder().setTransport(httpTransport)
        .setJsonFactory(jsonFactory)
        .setClientSecrets(clientId, clientSecret)
        .build()
        .setRefreshToken(refreshToken)
    )
  }
}

trait ClientSecrets {
  val clientId: String
  val clientSecret: String
}

final case class SimpleClientSecrets(clientId: String, clientSecret: String) extends ClientSecrets
