package drs.localizer.accesstokens

import cats.syntax.validated._
import com.google.auth.oauth2.GoogleCredentials
import common.validation.ErrorOr.ErrorOr

import scala.jdk.CollectionConverters._
import scala.util.{Failure, Success, Try}

/**
 * Strategy for obtaining an access token from Google Application Default credentials that are assumed to already exist.
 */
case object GoogleAccessTokenStrategy extends AccessTokenStrategy {
  private final val UserInfoEmailScope = "https://www.googleapis.com/auth/userinfo.email"
  private final val UserInfoProfileScope = "https://www.googleapis.com/auth/userinfo.profile"
  private final val UserInfoScopes = List(UserInfoEmailScope, UserInfoProfileScope)

  override def getAccessToken(): ErrorOr[String] = {
    Try {
      val scopedCredentials = GoogleCredentials.getApplicationDefault().createScoped(UserInfoScopes.asJava)
      scopedCredentials.refreshAccessToken().getTokenValue
    } match {
      case Success(null) => "null token value attempting to refresh access token".invalidNel
      case Success(value) => value.validNel
      case Failure(e) => s"Failed to refresh access token: ${e.getMessage}".invalidNel
    }
  }
}
