package drs.localizer

import cloud.nio.impl.drs.{DrsConfig, DrsPathResolver}
import com.google.auth.oauth2.GoogleCredentials
import org.apache.http.impl.client.HttpClientBuilder

import scala.collection.JavaConverters._

class LocalizerDrsPathResolver(drsConfig: DrsConfig,
                               httpClientBuilder: HttpClientBuilder) extends DrsPathResolver(drsConfig, httpClientBuilder) {

  private val UserInfoEmailScope = "https://www.googleapis.com/auth/userinfo.email"
  private val UserInfoProfileScope = "https://www.googleapis.com/auth/userinfo.profile"
  private val UserInfoScopes = List(UserInfoEmailScope, UserInfoProfileScope)

  override def getAccessToken: String = {
    val scopedCredentials = GoogleCredentials.getApplicationDefault().createScoped(UserInfoScopes.asJava)
    scopedCredentials.refreshAccessToken().getTokenValue
  }
}
