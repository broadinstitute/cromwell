package cromwell.util.google

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
import cromwell.engine.io.gcs._
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

object GoogleCredentialFactory extends GoogleCredentialFactory {
  // gcloudConf is a Try, calling .get will throw an exception containing all validation failures.
  override lazy val GoogleConf = GoogleConfiguration.gcloudConf.get

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

}

/**
  * Should not be instantiated except in tests. Use the GoogleCredentialFactory object instead.
  */
abstract class GoogleCredentialFactory {

  // Abstract to allow use of custom configuration in the tests.
  val GoogleConf: GoogleConfiguration

  private lazy val log = LoggerFactory.getLogger("GoogleCredentialFactory")
  lazy val jsonFactory = JacksonFactory.getDefaultInstance
  lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport

  lazy val fromCromwellAuthScheme: Credential = GoogleConf.cromwellAuthMode match {
    case user: UserMode => validateCredentials(forUser(user))
    case service: ServiceAccountMode => validateCredentials(forServiceAccount(service))
    case ApplicationDefaultMode => validateCredentials(forApplicationDefaultCredentials())
  }

  lazy val fromUserAuthScheme: (String) => Try[Credential] = (forClientSecrets _).andThen(_ map validateCredentials)

  private def validateCredentials(credential: Credential) = {
    Try(credential.refreshToken()) match {
      case Failure(ex) => throw new Throwable(s"Google credentials are invalid: ${ex.getMessage}")
      case Success(_) => credential
    }
  }

  private def filePathToSecrets(secrets: String, jsonFactory: JsonFactory) = {
    val secretsPath = Paths.get(secrets)
    if(!secretsPath.toFile.canRead) {
      log.warn("Secrets file does not exist or is not readable.")
    }
    val secretStream = new InputStreamReader(new FileInputStream(secretsPath.toFile))

    GoogleClientSecrets.load(jsonFactory, secretStream)
  }

  private def forUser(userConf: UserMode): Credential = {
    val user = userConf.user
    val clientSecrets = filePathToSecrets(userConf.secretsFile, jsonFactory)
    val dataStore = Paths.get(userConf.datastoreDir)
    val dataStoreFactory = new FileDataStoreFactory(dataStore.toFile)
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport,
      jsonFactory,
      clientSecrets,
      GoogleScopes.Scopes.asJava).setDataStoreFactory(dataStoreFactory).build
    new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(user)
  }

  private def forServiceAccount(serviceAccountMode: ServiceAccountMode): Credential = {
    val pemFile = new File(serviceAccountMode.pemPath)
    if (!pemFile.exists()) {
      throw new FileNotFoundException(s"Pem file ${Paths.get(serviceAccountMode.pemPath).toAbsolutePath} does not exist")
    }
    new GoogleCredential.Builder().setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setServiceAccountId(serviceAccountMode.accountId)
      .setServiceAccountScopes(GoogleScopes.Scopes.asJava)
      .setServiceAccountPrivateKeyFromPemFile(pemFile)
      .build()
  }

  private def forClientSecrets(token: String): Try[Credential] = {
    GoogleConf.userAuthMode collect {
      case refresh: Refresh =>
        Success(new GoogleCredential.Builder().setTransport(httpTransport)
          .setJsonFactory(jsonFactory)
          .setClientSecrets(refresh.clientSecrets.clientId, refresh.clientSecrets.clientSecret)
          .build()
          .setRefreshToken(token))
      case unrecognized => Failure(new IllegalArgumentException(s"Unrecognized userSchemeAuthentication: $unrecognized"))
    } getOrElse Failure(new Throwable("No user authentication configuration has been found."))
  }

  private def forApplicationDefaultCredentials(): Credential = {
    try {
      GoogleCredential.getApplicationDefault().createScoped(GoogleScopes.Scopes.asJava)
    } catch {
      case e: IOException =>
       log.warn("Failed to get application default credentials", e)
       throw e
    }
  }
}
