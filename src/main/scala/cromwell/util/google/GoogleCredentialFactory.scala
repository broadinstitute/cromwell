package cromwell.util.google

import java.io.{File, FileInputStream, InputStreamReader}
import java.nio.file.Paths

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp
import com.google.api.client.googleapis.auth.oauth2.{GoogleAuthorizationCodeFlow, GoogleClientSecrets, GoogleCredential}
import com.google.api.client.googleapis.extensions.java6.auth.oauth2.GooglePromptReceiver
import com.google.api.client.http.HttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.client.util.store.FileDataStoreFactory
import com.typesafe.config.{Config, ConfigFactory}
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._

object GoogleCredentialFactory {
  private lazy val GoogleConf = ConfigFactory.load.getConfig("google")
  private lazy val GoogleAuthScheme = GoogleConf.getString("authScheme").toLowerCase
  private lazy val log = LoggerFactory.getLogger("GoogleCredentialFactory")

  lazy val from: (JsonFactory, HttpTransport) => Credential = GoogleAuthScheme match {
    case "user" => forUser(GoogleConf.getConfig("userAuth"))
    case "service" => forServiceAccount(GoogleConf.getConfig("serviceAuth"))
  }

  private def filePathToSecrets(secrets: String, jsonFactory: JsonFactory) = {
    val secretsPath = Paths.get(secrets)
    if(!secretsPath.toFile.canRead) {
      log.warn("Secrets file does not exist or is not readable.")
    }
    val secretStream = new InputStreamReader(new FileInputStream(secretsPath.toFile))

    GoogleClientSecrets.load(jsonFactory, secretStream)
  }

  private def forUser(config: Config)(jsonFactory: JsonFactory, httpTransport: HttpTransport): Credential = {
    val user = config.getString("user")
    val clientSecrets = filePathToSecrets(config.getString("secretsFile"), jsonFactory)
    val dataStore = Paths.get(config.getString("dataStoreDir"))
    val dataStoreFactory = new FileDataStoreFactory(dataStore.toFile)
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport,
      jsonFactory,
      clientSecrets,
      GoogleScopes.Scopes.asJava).setDataStoreFactory(dataStoreFactory).build
    new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(user)
  }

  private def forServiceAccount(config: Config)(jsonFactory: JsonFactory, httpTransport: HttpTransport): Credential = {
    new GoogleCredential.Builder().setTransport(httpTransport)
      .setJsonFactory(jsonFactory)
      .setServiceAccountId(config.getString("serviceAccountId"))
      .setServiceAccountScopes(GoogleScopes.Scopes.asJava)
      .setServiceAccountPrivateKeyFromP12File(new File(config.getString("p12File")))
    //  .setServiceAccountUser(GoogleUser) FIXME: Dig into how impersonation works and if we even care
      .build()
  }
}
