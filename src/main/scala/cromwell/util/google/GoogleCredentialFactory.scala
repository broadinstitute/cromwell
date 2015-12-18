package cromwell.util.google

import java.io.FileReader
import java.nio.file.{Path, Paths}

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp
import com.google.api.client.googleapis.auth.oauth2.{GoogleAuthorizationCodeFlow, GoogleClientSecrets, GoogleCredential}
import com.google.api.client.googleapis.extensions.java6.auth.oauth2.GooglePromptReceiver
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.client.util.store.FileDataStoreFactory
import com.typesafe.config.{Config, ConfigFactory}

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

/**
  * Holds a reference to a Google Credential, optionally refreshing the access token before returning the Credential.
  */
trait GoogleCredentialFactory {
  lazy val rawCredential = initCredential()

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
  def freshCredential: Try[Credential] = {
    val stillValid = Option(rawCredential.getExpiresInSeconds).exists(_ > 60)
    if (stillValid || rawCredential.refreshToken()) {
      Success(rawCredential)
    } else {
      Failure(new Exception("Unable to refresh token"))
    }
  }

  /** Returns the initial credential that will be refreshed. */
  protected def initCredential(): Credential

  protected lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport
  protected lazy val jsonFactory = JacksonFactory.getDefaultInstance
  protected lazy val scopes = GoogleScopes.Scopes.asJava
}

object GoogleCredentialFactory {
  def fromAuthScheme(config: Config) = {
    config.getString("google.authScheme").toLowerCase match {
      case "user" => new UserCredentialFactory(config)
      case "service" => new ServiceAccountCredentialFactory(config)
    }
  }

  lazy val fromAuthScheme: Credential = validateCredentials(
    GoogleCredentialFactory.fromAuthScheme(ConfigFactory.load).rawCredential)

  lazy val forRefreshToken: (ClientSecrets, String) => Credential = forClientSecrets

  private def forClientSecrets(secrets: ClientSecrets, token: String): Credential = {
    validateCredentials(new RefreshTokenCredentialFactory(secrets.clientId, secrets.clientSecret, token).rawCredential)
  }

  private def validateCredentials(credential: Credential) = {
    Try(credential.refreshToken()) match {
      case Failure(ex) => throw new Throwable(s"Google credentials are invalid: ${ex.getMessage}")
      case Success(_) => credential
    }
  }
}

/**
  * Creates credentials for a userId.
  * @param userId The user id to authenticate as.
  * @param secretsPath The path to a file that holds the user's secrets.
  * @param dataStoreDir A directory that holds the possibly existing cached credentials for the user's credentials.
  */
class UserCredentialFactory(userId: String, secretsPath: Path, dataStoreDir: Path) extends GoogleCredentialFactory {
  def this(config: Config) = {
    this(
      config.getString("google.userAuth.user"),
      Paths.get(config.getString("google.userAuth.secretsFile")),
      Paths.get(config.getString("google.userAuth.dataStoreDir"))
    )
  }

  protected def initCredential(): Credential = {
    val clientSecrets = GoogleClientSecrets.load(jsonFactory, new FileReader(secretsPath.toFile))
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport, jsonFactory, clientSecrets, scopes)
      .setDataStoreFactory(new FileDataStoreFactory(dataStoreDir.toFile))
      .build()
    new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(userId)
  }
}

/**
  * Creates a credential for a service account.
  * @param serviceAccountId The id of the service account.
  * @param pemPath The path to the pem file for the service account.
  */
class ServiceAccountCredentialFactory(serviceAccountId: String, pemPath: Path) extends GoogleCredentialFactory {
  def this(config: Config) = {
    this(
      config.getString("google.serviceAuth.serviceAccountId"),
      Paths.get(config.getString("google.serviceAuth.pemFile"))
    )
  }

  protected def initCredential(): Credential = {
    new GoogleCredential.Builder()
      .setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setServiceAccountScopes(scopes)
      .setServiceAccountId(serviceAccountId)
      .setServiceAccountPrivateKeyFromPemFile(pemPath.toFile)
      .build()
  }
}

/**
  * Creates a credential for a refresh token.
  * @param clientId The client id.
  * @param clientSecret The client secret.
  * @param refreshToken The refresh token.
  */
class RefreshTokenCredentialFactory(clientId: String, clientSecret: String,
                                    refreshToken: String) extends GoogleCredentialFactory {
  def this(config: Config, refreshToken: String) = {
    this(
      config.getString("google.localizeWithRefreshToken.client_id"),
      config.getString("google.localizeWithRefreshToken.client_secret"),
      refreshToken
    )
  }

  protected def initCredential(): Credential = {
    new GoogleCredential.Builder()
      .setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setClientSecrets(clientId, clientSecret)
      .build()
      .setRefreshToken(refreshToken)
  }
}
