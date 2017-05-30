package cromwell.filesystems.gcs.auth

import java.io.FileNotFoundException

import better.files._
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.storage.StorageScopes
import com.google.auth.Credentials
import com.google.auth.http.HttpTransportFactory
import com.google.auth.oauth2.{GoogleCredentials, ServiceAccountCredentials, UserCredentials}
import com.google.cloud.NoCredentials
import cromwell.core.WorkflowOptions
import cromwell.filesystems.gcs.auth.GoogleAuthMode._
import cromwell.filesystems.gcs.auth.ServiceAccountMode.{CredentialFileFormat, JsonFileFormat, PemFileFormat}
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

object GoogleAuthMode {

  lazy val jsonFactory = JacksonFactory.getDefaultInstance
  lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport
  lazy val HttpTransportFactory = new HttpTransportFactory {
    override def create() = {
      httpTransport
    }
  }

  val RefreshTokenOptionKey = "refresh_token"
  val GcsScopes = List(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE
  ).asJava

  def checkReadable(file: File) = {
    if (!file.isReadable) throw new FileNotFoundException(s"File $file does not exist or is not readable")
  }

  case object MockAuthMode extends GoogleAuthMode {
    override def name = "no_auth"
    override def credential(options: WorkflowOptions): Credentials = NoCredentials.getInstance()
  }
}


sealed trait GoogleAuthMode {
  protected lazy val log = LoggerFactory.getLogger(getClass.getSimpleName)

  /**
    * Validate the auth mode against provided options
    */
  def validate(options: WorkflowOptions): Unit = {()}

  def name: String
  // Create a Credential object from the google.api.client.auth library (https://github.com/google/google-api-java-client)
  def credential(options: WorkflowOptions): Credentials

  def requiresAuthFile: Boolean = false

  protected def validateCredential(credential: Credentials) = validate(credential, () => credential.refresh())

  private def validate[T](credential: T, validation: () => Any): T = {
    Try(validation()) match {
      case Failure(ex) => throw new RuntimeException(s"Google credentials are invalid: ${ex.getMessage}", ex)
      case Success(_) => credential
    }
  }
}

object ServiceAccountMode {
  sealed trait CredentialFileFormat {
    def file: String
  }
  case class PemFileFormat(accountId: String, file: String) extends CredentialFileFormat
  case class JsonFileFormat(file: String) extends CredentialFileFormat
}

final case class ServiceAccountMode(override val name: String,
                                    fileFormat: CredentialFileFormat,
                                    scopes: java.util.List[String]) extends GoogleAuthMode {
  private val credentialsFile = File(fileFormat.file)
  checkReadable(credentialsFile)

  private lazy val _credential: Credentials = {
    val serviceAccount = fileFormat match {
      case PemFileFormat(accountId, _) => 
        log.warn("The PEM file format will be deprecated in the upcoming Cromwell version. Please use JSON instead.")
        ServiceAccountCredentials.fromPkcs8(accountId, accountId, credentialsFile.contentAsString, null, scopes)
      case _: JsonFileFormat => ServiceAccountCredentials.fromStream(credentialsFile.newInputStream).createScoped(scopes)
    }
    
    validateCredential(serviceAccount)
  }

  override def credential(options: WorkflowOptions): Credentials = _credential
}

final case class UserMode(override val name: String,
                          user: String,
                          secretsPath: String,
                          datastoreDir: String,
                          scopes: java.util.List[String]) extends GoogleAuthMode {

  private lazy val secretsStream = {
    val secretsFile = File(secretsPath)
    checkReadable(secretsFile)
    secretsFile.newInputStream
  }

  private lazy val _credential: Credentials = {
    validateCredential(UserCredentials.fromStream(secretsStream))
  }

  override def credential(options: WorkflowOptions) = _credential
}

private object ApplicationDefault {
  private [auth] lazy val _Credential: Credentials = GoogleCredentials.getApplicationDefault
}

final case class ApplicationDefaultMode(name: String) extends GoogleAuthMode {
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

  override def credential(options: WorkflowOptions): Credentials = {
    val refreshToken = extractRefreshToken(options)
    validateCredential(
      new UserCredentials(clientId, clientSecret, refreshToken, null, GoogleAuthMode.HttpTransportFactory, null)
    )
  }
}

trait ClientSecrets {
  val clientId: String
  val clientSecret: String
}

final case class SimpleClientSecrets(clientId: String, clientSecret: String) extends ClientSecrets
