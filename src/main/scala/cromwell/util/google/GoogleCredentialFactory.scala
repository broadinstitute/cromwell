package cromwell.util.google

import java.io.{File, FileInputStream, InputStreamReader}
import java.nio.file.{Paths, Path}
import java.security.PrivateKey

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp
import com.google.api.client.googleapis.auth.oauth2.GoogleCredential.Builder
import com.google.api.client.googleapis.auth.oauth2.{GoogleCredential, GoogleAuthorizationCodeFlow, GoogleClientSecrets}
import com.google.api.client.googleapis.extensions.java6.auth.oauth2.GooglePromptReceiver
import com.google.api.client.http.HttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.client.util.store.FileDataStoreFactory
import com.google.api.services.storage.StorageScopes
import com.typesafe.config.{Config, ConfigFactory}
import scala.collection.JavaConverters._

// FIXME: Use Tex's service account

object GoogleCredentialFactory {
  private lazy val GoogleConf = ConfigFactory.load.getConfig("google")
  private lazy val GoogleUser = GoogleConf.getString("user")
  private lazy val GoogleAuthScheme = GoogleConf.getString("authScheme").toLowerCase

  lazy val from: (JsonFactory, HttpTransport) => Credential = GoogleAuthScheme match {
    case "user" => forUser(GoogleConf.getConfig("userAuth"))
    case "service" => forServiceAccount(GoogleConf.getConfig("serviceAuth"))
  }

  private def forUser(config: Config)(jsonFactory: JsonFactory, httpTransport: HttpTransport): Credential = {
    val secrets = Paths.get(config.getString("secretsFile"))
    val secretStream = new InputStreamReader(new FileInputStream(secrets.toFile))
    val clientSecrets = GoogleClientSecrets.load(jsonFactory, secretStream)
    val dataStore = Paths.get(config.getString("dataStoreDir"))
    val dataStoreFactory = new FileDataStoreFactory(dataStore.toFile)
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport,
      jsonFactory,
      clientSecrets,
      GoogleScopes.Scopes.asJava).setDataStoreFactory(dataStoreFactory).build
    new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(GoogleUser)
  }

  private def forServiceAccount(config: Config)(jsonFactory: JsonFactory, httpTransport: HttpTransport): Credential = ???
}
