package cromwell.filesystems.gcs.auth

import com.google.api.client.auth.oauth2.Credential
import com.google.cloud.AuthCredentials

case class GoogleCredentialBundle(credential: Credential, authCredential: AuthCredentials)
