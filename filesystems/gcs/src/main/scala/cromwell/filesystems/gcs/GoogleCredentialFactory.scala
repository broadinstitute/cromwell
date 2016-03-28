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
import org.slf4j.LoggerFactory

import scala.util.{Failure, Success, Try}

/**
  * Should not be instantiated except in tests. Use the GoogleCredentialFactory object instead.
  */
object GoogleCredentialFactory {

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

  def apply(conf: GoogleAuthMode, scopes: GoogleScopes, refreshToken: Option[RefreshToken] = None) = conf match {
    case user: UserMode => validateCredentials(forUser(user, scopes))
    case service: ServiceAccountMode => validateCredentials(forServiceAccount(service, scopes))
    case refresh: RefreshTokenMode => validateCredentials(refreshTokenMode(refresh, refreshToken, scopes))
    case ApplicationDefaultMode => validateCredentials(forApplicationDefaultCredentials(scopes))
  }

  private lazy val log = LoggerFactory.getLogger("GoogleCredentialFactory")
  lazy val jsonFactory = JacksonFactory.getDefaultInstance
  lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport

  private def refreshTokenMode(refreshTokenMode: RefreshTokenMode, refreshToken: Option[RefreshToken], scopes: GoogleScopes) = {
    refreshToken map forClientSecrets(refreshTokenMode, scopes) getOrElse {
      throw new IllegalArgumentException(
        "Could not build Google Credentials: A refresh token needs to be provided when using refresh token authentication mode."
      )
    }
  }

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

  private def forUser(userConf: UserMode, scopes: GoogleScopes): Credential = {
    val user = userConf.user
    val clientSecrets = filePathToSecrets(userConf.secretsFile, jsonFactory)
    val dataStore = Paths.get(userConf.datastoreDir)
    val dataStoreFactory = new FileDataStoreFactory(dataStore.toFile)
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport,
      jsonFactory,
      clientSecrets,
      scopes).setDataStoreFactory(dataStoreFactory).build
    new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(user)
  }

  private def forServiceAccount(serviceAccountMode: ServiceAccountMode, scopes: GoogleScopes): Credential = {
    val pemFile = new File(serviceAccountMode.pemPath)
    if (!pemFile.exists()) {
      throw new FileNotFoundException(s"Pem file ${Paths.get(serviceAccountMode.pemPath).toAbsolutePath} does not exist")
    }
    new GoogleCredential.Builder().setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setServiceAccountId(serviceAccountMode.accountId)
      .setServiceAccountScopes(scopes)
      .setServiceAccountPrivateKeyFromPemFile(pemFile)
      .build()
  }

  private def forClientSecrets(conf: RefreshTokenMode, scopes: GoogleScopes)(token: String): Credential = {
    new GoogleCredential.Builder().setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setClientSecrets(conf.clientSecrets.clientId, conf.clientSecrets.clientSecret)
      .build()
      .createScoped(scopes)
      .setRefreshToken(token)
  }

  private def forApplicationDefaultCredentials(scopes: GoogleScopes): Credential = {
    try {
      GoogleCredential.getApplicationDefault().createScoped(scopes)
    } catch {
      case e: IOException =>
       log.warn("Failed to get application default credentials", e)
       throw e
    }
  }
}
