package cromwell.filesystems.gcs.auth

import java.io.Serializable
import java.util.Objects

import com.google.auth.oauth2.{ClientId, GoogleCredentials, UserCredentials}
import com.google.cloud.{AuthCredentials, RestorableState}

class RefreshableOAuth2Credentials(refreshToken: String, clientId: ClientId) extends AuthCredentials {
  private val _credentials: GoogleCredentials = new UserCredentials(clientId.getClientId, clientId.getClientSecret, refreshToken)

  private class RefreshableOAuth2CredentialsState(val refreshToken: String, val clientId: ClientId) extends RestorableState[AuthCredentials] with Serializable {

    override def restore: AuthCredentials = new RefreshableOAuth2Credentials(refreshToken, clientId)

    override def hashCode: Int = Objects.hash(refreshToken, clientId.getClientId, clientId.getClientSecret)

    override def equals(obj: Any): Boolean = {
      obj.isInstanceOf[RefreshableOAuth2CredentialsState] && {
        val other = obj.asInstanceOf[RefreshableOAuth2CredentialsState]
        Objects.equals(refreshToken, other.refreshToken) &&
        Objects.equals(clientId.getClientId, other.clientId.getClientId) &&
        Objects.equals(clientId.getClientSecret, other.clientId.getClientSecret)
      }
    }
  }

  override def credentials: GoogleCredentials = _credentials

  def capture: RestorableState[AuthCredentials] = new RefreshableOAuth2CredentialsState(refreshToken, clientId)
}
