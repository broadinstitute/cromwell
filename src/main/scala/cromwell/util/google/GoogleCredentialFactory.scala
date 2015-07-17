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
import com.typesafe.config.ConfigFactory
import scala.collection.JavaConverters._

// FIXME: Use Tex's service account

object GoogleCredentialFactory {
  private lazy val GoogleConf = ConfigFactory.load.getConfig("google")
  lazy val GoogleSecrets = Paths.get(GoogleConf.getString("secretsFile"))
  lazy val GoogleUser = GoogleConf.getString("user")

  def from(jsonFactory: JsonFactory, httpTransport: HttpTransport): Credential = {
    val secretStream = new InputStreamReader(new FileInputStream(GoogleSecrets.toFile))
    val clientSecrets = GoogleClientSecrets.load(jsonFactory, secretStream)
    // FIXME: The following shouldn't be hardcoded
    val dataStoreFactory = new FileDataStoreFactory(new File(System.getProperty("user.home"), ".jes-google-alpha"))
    val flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport,
      jsonFactory,
      clientSecrets,
      GoogleScopes.Scopes.asJava).setDataStoreFactory(dataStoreFactory).build
    new AuthorizationCodeInstalledApp(flow, new GooglePromptReceiver).authorize(GoogleUser)
  }
}
