package cromwell.cloudsupport.gcp.auth

import java.io.{ByteArrayInputStream, FileNotFoundException}
import java.net.HttpURLConnection._

import better.files.File
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.http.HttpResponseException
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.storage.StorageScopes
import com.google.auth.Credentials
import com.google.auth.http.HttpTransportFactory
import com.google.auth.oauth2.{GoogleCredentials, ServiceAccountCredentials, UserCredentials}
import com.google.cloud.NoCredentials
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode._
import cromwell.cloudsupport.gcp.auth.ServiceAccountMode.{CredentialFileFormat, JsonFileFormat, PemFileFormat}
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

object GoogleAuthMode {

  type OptionLookup = String => String
  lazy val jsonFactory = JacksonFactory.getDefaultInstance
  lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport
  lazy val HttpTransportFactory = new HttpTransportFactory {
    override def create() = {
      httpTransport
    }
  }

  val RefreshTokenOptionKey = "refresh_token"
  val UserServiceAccountKey = "user_service_account_json"

  val GcsScopes = List(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE
  ).asJava

  def checkReadable(file: File) = {
    if (!file.isReadable) throw new FileNotFoundException(s"File $file does not exist or is not readable")
  }

  def isFatal(ex: Throwable) = {
    ex match {
      case http: HttpResponseException =>
        // Using HttpURLConnection fields as com.google.api.client.http.HttpStatusCodes doesn't have Bad Request (400)
        http.getStatusCode == HTTP_UNAUTHORIZED ||
          http.getStatusCode == HTTP_FORBIDDEN ||
          http.getStatusCode == HTTP_BAD_REQUEST
      case _: OptionLookupException => true
      case _ => false
    }
  }

  def extract(options: OptionLookup, key: String): String = {
    Try(options(key)) match {
      case Success(result) => result
      case Failure(throwable) => throw new OptionLookupException(key, throwable)
    }
  }

}


sealed trait GoogleAuthMode {
  protected lazy val log = LoggerFactory.getLogger(getClass.getSimpleName)

  /**
    * Validate the auth mode against provided options
    */
  def validate(options: OptionLookup): Unit = {
    ()
  }

  def name: String

  // Create a Credential object from the google.api.client.auth library (https://github.com/google/google-api-java-client)
  def credential(options: OptionLookup): Credentials

  def requiresAuthFile: Boolean = false

  /**
    * Enables swapping out credential validation for various testing purposes ONLY.
    *
    * All traits in this file are sealed, all classes final, meaning things like Mockito or other java/scala overrides
    * cannot work.
    */
  private[auth] var credentialValidation: (Credentials => Unit) = credentials => credentials.refresh()

  protected def validateCredential(credential: Credentials) = {
    Try(credentialValidation(credential)) match {
      case Failure(ex) => throw new RuntimeException(s"Google credentials are invalid: ${ex.getMessage}", ex)
      case Success(_) => credential
    }
  }
}

case object MockAuthMode extends GoogleAuthMode {
  override val name = "no_auth"

  override def credential(options: OptionLookup): Credentials = NoCredentials.getInstance
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

    // Validate credentials synchronously here, without retry.
    // It's very unlikely to fail as it should not happen more than a few times
    // (one for the engine and for each backend using it) per Cromwell instance.
    validateCredential(serviceAccount)
  }

  override def credential(options: OptionLookup): Credentials = _credential
}

final case class UserServiceAccountMode(override val name: String, scopes: java.util.List[String]) extends GoogleAuthMode {
  private def extractServiceAccount(options: OptionLookup): String = {
    extract(options, UserServiceAccountKey)
  }

  override def validate(options: OptionLookup): Unit = {
    extractServiceAccount(options)
    ()
  }

  override def credential(options: OptionLookup): Credentials = {
    ServiceAccountCredentials
      .fromStream(new ByteArrayInputStream(extractServiceAccount(options).getBytes("UTF-8")))
      .createScoped(scopes)
  }
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

  override def credential(options: OptionLookup): Credentials = _credential
}

object ApplicationDefaultMode {
  private lazy val _credential: Credentials = GoogleCredentials.getApplicationDefault
}

final case class ApplicationDefaultMode(name: String) extends GoogleAuthMode {
  override def credential(options: OptionLookup): Credentials = ApplicationDefaultMode._credential
}

final case class RefreshTokenMode(name: String,
                                  clientId: String,
                                  clientSecret: String,
                                  scopes: java.util.List[String]) extends GoogleAuthMode with ClientSecrets {

  import GoogleAuthMode._

  override def requiresAuthFile = true

  private def extractRefreshToken(options: OptionLookup): String = {
    extract(options, RefreshTokenOptionKey)
  }

  override def validate(options: OptionLookup) = {
    extractRefreshToken(options)
    ()
  }

  override def credential(options: OptionLookup): Credentials = {
    val refreshToken = extractRefreshToken(options)
    validateCredential(
      UserCredentials
        .newBuilder()
        .setClientId(clientId)
        .setClientSecret(clientSecret)
        .setRefreshToken(refreshToken)
        .setHttpTransportFactory(GoogleAuthMode.HttpTransportFactory)
        .build()
    )
  }
}

sealed trait ClientSecrets {
  val clientId: String
  val clientSecret: String
}

final case class SimpleClientSecrets(clientId: String, clientSecret: String) extends ClientSecrets

class OptionLookupException(val key: String, cause: Throwable) extends RuntimeException(key, cause)
