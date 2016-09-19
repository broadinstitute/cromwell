package cromwell.backend.impl.jes.authentication

import com.google.api.client.auth.oauth2.Credential

case class JesCredentials(genomicsCredential: Credential, gcsCredential: Credential)
