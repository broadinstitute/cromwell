package cromwell.cloudsupport.gcp.auth

import java.io.{ByteArrayInputStream, FileNotFoundException, InputStream}
import java.net.HttpURLConnection._

import better.files.File
import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.googleapis.testing.auth.oauth2.MockGoogleCredential
import com.google.api.client.http.HttpResponseException
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.cloudkms.v1.CloudKMS
import com.google.api.services.compute.ComputeScopes
import com.google.api.services.genomics.v2alpha1.GenomicsScopes
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
  val NoOptionLookup = noOptionLookup _

  /**
    * Enables swapping out credential validation for various testing purposes ONLY.
    *
    * All traits in this file are sealed, all classes final, meaning things like Mockito or other java/scala overrides
    * cannot work.
    */
  type CredentialsValidation = Credentials => Unit
  private[auth] val NoCredentialsValidation = mouse.ignore _

  private def noOptionLookup(string: String): Nothing = {
    throw new UnsupportedOperationException(s"cannot lookup $string")
  }

  lazy val jsonFactory = JacksonFactory.getDefaultInstance
  lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport
  lazy val HttpTransportFactory = new HttpTransportFactory {
    override def create() = {
      httpTransport
    }
  }

  val RefreshTokenOptionKey = "refresh_token"
  val UserServiceAccountKey = "user_service_account_json"
  val DockerCredentialsEncryptionKeyNameKey = "docker_credentials_key_name"
  val DockerCredentialsTokenKey = "docker_credentials_token"

  private val PipelinesApiScopes = List(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE,
    GenomicsScopes.GENOMICS,
    ComputeScopes.COMPUTE
  )

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

  def encryptKms(keyName: String, credential: GoogleCredential, plainText: String) = {
    import com.google.api.services.cloudkms.v1.CloudKMSScopes

    // Depending on the environment that provides the default credentials (e.g. Compute Engine, App
    // Engine), the credentials may require us to specify the scopes we need explicitly.
    // Check for this case, and inject the scope if required.
    val scopedCredential = if (credential.createScopedRequired) credential.createScoped(CloudKMSScopes.all) else credential

    val kms = new CloudKMS.Builder(httpTransport, jsonFactory, scopedCredential)
      .setApplicationName("cromwell")
      .build()

    import com.google.api.services.cloudkms.v1.model.EncryptRequest
    val request = new EncryptRequest().encodePlaintext(plainText.toCharArray.map(_.toByte))
    val response = kms.projects.locations.keyRings.cryptoKeys.encrypt(keyName, request).execute
    response.getCiphertext
  }

  /** Used for both checking that the credential is valid and creating a fresh credential. */
  private def refreshCredentials(credentials: Credentials): Unit = {
    credentials.refresh()
  }
}


sealed trait GoogleAuthMode {
  protected lazy val log = LoggerFactory.getLogger(getClass.getSimpleName)

  def name: String

  // Create a Credential object from the google.api.client.auth library (https://github.com/google/google-api-java-client)
  private[auth] def credentials(options: OptionLookup, scopes: java.util.Collection[String]): Credentials

  /**
    * Create a credential object suitable for use with Pipelines API.
    *
    * @param options A lookup for external credential information.
    * @return Credentials with scopes compatible with the Genomics API compute and storage.
    */
  def pipelinesApiCredentials(options: OptionLookup): Credentials = {
    credentials(options, PipelinesApiScopes.asJavaCollection)
  }

  /**
    * Alias for credentials(GoogleAuthMode.NoOptionLookup, scopes).
    * Only valid for credentials that are NOT externally provided, such as ApplicationDefault.
    */
  def credentials(scopes: Iterable[String]): Credentials = {
    credentials(GoogleAuthMode.NoOptionLookup, scopes.asJavaCollection)
  }

  /**
    * Alias for credentials(GoogleAuthMode.NoOptionLookup, scopes).
    * Only valid for credentials that are NOT externally provided, such as ApplicationDefault.
    */
  def credentials(scopes: java.util.Collection[String]): Credentials = {
    credentials(GoogleAuthMode.NoOptionLookup, scopes)
  }

  /**
    * Alias for credentials(GoogleAuthMode.NoOptionLookup, Set.empty).
    * Only valid for credentials that are NOT externally provided and do not need scopes, such as ApplicationDefault.
    */
  private[auth] def credentials(): Credentials = {
    credentials(GoogleAuthMode.NoOptionLookup, java.util.Collections.emptySet[String])
  }

  /**
    * Alias for credentials(options, Set.empty).
    * Only valid for credentials that are NOT externally provided and do not need scopes, such as ApplicationDefault.
    */
  private[auth] def credentials(options: OptionLookup): Credentials = {
    credentials(options, java.util.Collections.emptySet[String])
  }

  def requiresAuthFile: Boolean = false

  /**
    * Enables swapping out credential validation for various testing purposes ONLY.
    *
    * All traits in this file are sealed, all classes final, meaning things like Mockito or other java/scala overrides
    * cannot work.
    */
  private[auth] var credentialsValidation: CredentialsValidation = refreshCredentials

  protected def validateCredentials[A <: Credentials](credential: A): credential.type = {
    Try(credentialsValidation(credential)) match {
      case Failure(ex) => throw new RuntimeException(s"Google credentials are invalid: ${ex.getMessage}", ex)
      case Success(_) => credential
    }
  }

  def apiClientGoogleCredential(options: OptionLookup): Option[GoogleCredential] = None
}

case object MockAuthMode extends GoogleAuthMode {
  override val name = "no_auth"

  override def credentials(unusedOptions: OptionLookup, unusedScopes: java.util.Collection[String]): NoCredentials = {
    NoCredentials.getInstance
  }

  override def apiClientGoogleCredential(options: OptionLookup): Option[MockGoogleCredential] = {
    Option(new MockGoogleCredential.Builder().build())
  }
}

object ServiceAccountMode {

  sealed trait CredentialFileFormat {
    def file: String
  }

  case class PemFileFormat(accountId: String, file: String) extends CredentialFileFormat

  case class JsonFileFormat(file: String) extends CredentialFileFormat

}

trait HasApiClientGoogleCredentialStream { self: GoogleAuthMode =>
  protected def credentialStream(options: OptionLookup): InputStream

  override def apiClientGoogleCredential(options: OptionLookup): Option[GoogleCredential] = Option(GoogleCredential.fromStream(credentialStream(options)))
}

final case class ServiceAccountMode(override val name: String,
                                    fileFormat: CredentialFileFormat)
  extends GoogleAuthMode with HasApiClientGoogleCredentialStream {
  private val credentialsFile = File(fileFormat.file)
  checkReadable(credentialsFile)

  override protected def credentialStream(options: OptionLookup): InputStream = credentialsFile.newInputStream

  private lazy val serviceAccountCredentials: ServiceAccountCredentials = {
    fileFormat match {
      case PemFileFormat(accountId, _) =>
        log.warn("The PEM file format will be deprecated in the upcoming Cromwell version. Please use JSON instead.")
        ServiceAccountCredentials.fromPkcs8(accountId, accountId, credentialsFile.contentAsString, null, null)
      case _: JsonFileFormat => ServiceAccountCredentials.fromStream(credentialsFile.newInputStream)
    }
  }

  override def credentials(unusedOptions: OptionLookup,
                           scopes: java.util.Collection[String]): GoogleCredentials = {
    val scopedCredentials = serviceAccountCredentials.createScoped(scopes)
    validateCredentials(scopedCredentials)
  }
}

final case class UserServiceAccountMode(override val name: String)
  extends GoogleAuthMode with HasApiClientGoogleCredentialStream {
  private def extractServiceAccount(options: OptionLookup): String = {
    extract(options, UserServiceAccountKey)
  }

  override protected def credentialStream(options: OptionLookup): InputStream = {
    new ByteArrayInputStream(extractServiceAccount(options).getBytes("UTF-8"))
  }

  override def credentials(options: OptionLookup, scopes: java.util.Collection[String]): GoogleCredentials = {
    val newCredentials = ServiceAccountCredentials.fromStream(credentialStream(options))
    val scopedCredentials: GoogleCredentials = newCredentials.createScoped(scopes)
    validateCredentials(scopedCredentials)
  }
}


final case class UserMode(override val name: String,
                          user: String,
                          secretsPath: String,
                          datastoreDir: String) extends GoogleAuthMode {

  private lazy val secretsStream = {
    val secretsFile = File(secretsPath)
    checkReadable(secretsFile)
    secretsFile.newInputStream
  }

  private lazy val userCredentials: UserCredentials = {
    validateCredentials(UserCredentials.fromStream(secretsStream))
  }

  override def credentials(unusedOptions: OptionLookup, unusedScopes: java.util.Collection[String]): Credentials = {
    userCredentials
  }
}

object ApplicationDefaultMode {
  private lazy val applicationDefaultCredentials: GoogleCredentials = GoogleCredentials.getApplicationDefault
}

final case class ApplicationDefaultMode(name: String) extends GoogleAuthMode {
  override def credentials(unusedOptions: OptionLookup,
                           unusedScopes: java.util.Collection[String]): GoogleCredentials = {
    ApplicationDefaultMode.applicationDefaultCredentials
  }

  override def apiClientGoogleCredential(unused: OptionLookup): Option[GoogleCredential] = Option(GoogleCredential.getApplicationDefault(httpTransport, jsonFactory))
}

final case class RefreshTokenMode(name: String,
                                  clientId: String,
                                  clientSecret: String) extends GoogleAuthMode with ClientSecrets {

  import GoogleAuthMode._

  override def requiresAuthFile = true

  private def extractRefreshToken(options: OptionLookup): String = {
    extract(options, RefreshTokenOptionKey)
  }

  override def credentials(options: OptionLookup, unusedScopes: java.util.Collection[String]): UserCredentials = {
    val refreshToken = extractRefreshToken(options)
    val newCredentials: UserCredentials = UserCredentials
      .newBuilder()
      .setClientId(clientId)
      .setClientSecret(clientSecret)
      .setRefreshToken(refreshToken)
      .setHttpTransportFactory(GoogleAuthMode.HttpTransportFactory)
      .build()
    validateCredentials(newCredentials)
  }
}

sealed trait ClientSecrets {
  val clientId: String
  val clientSecret: String
}

final case class SimpleClientSecrets(clientId: String, clientSecret: String) extends ClientSecrets

class OptionLookupException(val key: String, cause: Throwable) extends RuntimeException(key, cause)
