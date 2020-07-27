package drs.localizer

import cloud.nio.impl.drs.{DrsConfig, DrsPathResolver}
import com.google.auth.oauth2.GoogleCredentials
import drs.localizer.DrsLocalizerMain.UserInfoScopes
import org.apache.http.impl.client.HttpClientBuilder

import scala.collection.JavaConverters._

class LocalizerDrsPathResolver(drsConfig: DrsConfig,
                               httpClientBuilder: HttpClientBuilder) extends DrsPathResolver(drsConfig, httpClientBuilder) {

  override def getAccessToken: String = {
    val scopedCredentials = GoogleCredentials.getApplicationDefault().createScoped(UserInfoScopes.asJava)
    scopedCredentials.refreshAccessToken().getTokenValue
  }
}
