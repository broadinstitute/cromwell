package cromwell.backend.impl.jes

import com.google.api.client.googleapis.testing.auth.oauth2.MockGoogleCredential
import cromwell.backend.impl.jes.authentication.JesCredentials

object MockObjects {
  val mockCredential = new MockGoogleCredential.Builder().build()
  val mockCredentials = JesCredentials(mockCredential, mockCredential)
}
