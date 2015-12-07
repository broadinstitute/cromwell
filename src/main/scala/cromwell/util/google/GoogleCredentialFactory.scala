package cromwell.util.google

import java.io.{File, FileInputStream, InputStreamReader}
import java.nio.file.Paths

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp
import com.google.api.client.googleapis.auth.oauth2.{GoogleAuthorizationCodeFlow, GoogleClientSecrets, GoogleCredential}
import com.google.api.client.googleapis.extensions.java6.auth.oauth2.GooglePromptReceiver
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.client.util.store.FileDataStoreFactory
import com.typesafe.config.{Config, ConfigFactory}
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

object GoogleCredentialFactory {
  private lazy val GoogleConf = ConfigFactory.load.getConfig("google")
  private lazy val GoogleAuthScheme = GoogleConf.getString("authScheme").toLowerCase
  private lazy val log = LoggerFactory.getLogger("GoogleCredentialFactory")
  lazy val jsonFactory = JacksonFactory.getDefaultInstance
  lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport

  lazy val fromAuthScheme: Credential = GoogleAuthScheme match {
    case "user" => forUser(GoogleConf.getConfig("userAuth"))
    case "service" => forServiceAccount(GoogleConf.getConfig("serviceAuth"))
  }

  lazy val forRefreshToken: (ClientSecrets, String) => Credential = forClientSecrets

  private def filePathToSecrets(secrets: String, jsonFactory: JsonFactory) = {
    val secretsPath = Paths.get(secrets)
    if(!secretsPath.toFile.canRead) {
      log.warn("Secrets file does not exist or is not readable.")
    }
    val secretStream = new InputStreamReader(new FileInputStream(secretsPath.toFile))

    GoogleClientSecrets.load(jsonFactory, secretStream)
  }

  private def validateCredentials(credential: Credential) = {
    Try(credential.refreshToken()) match {
      case Failure(ex) => throw new Throwable(s"Google credentials are invalid: ${ex.getMessage}")
      case Success(_) => credential
    }
  }

  private def forUser(config: Config): Credential = {
    val user = config.getString("user")
    val clientSecrets = filePathToSecrets(config.getString("secretsFile"), jsonFactory)
    val dataStore = Paths.get(config.getString("dataStoreDir"))
    val dataStoreFactory = new FileDataStoreFactory(dataStore.toFile)
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport,
      jsonFactory,
      clientSecrets,
      GoogleScopes.Scopes.asJava).setDataStoreFactory(dataStoreFactory).build
    validateCredentials(new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(user))
  }

  private def forServiceAccount(config: Config): Credential = {
    validateCredentials(new GoogleCredential.Builder().setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setServiceAccountId(config.getString("serviceAccountId"))
      .setServiceAccountScopes(GoogleScopes.Scopes.asJava)
      .setServiceAccountPrivateKeyFromPemFile(new File(config.getString("pemFile")))
      .build())
  }

  private def forClientSecrets(secrets: ClientSecrets, token: String): Credential = {
    validateCredentials(new GoogleCredential.Builder().setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setClientSecrets(secrets.clientId, secrets.clientSecret)
      .build()
      .setRefreshToken(token))
  }
}
