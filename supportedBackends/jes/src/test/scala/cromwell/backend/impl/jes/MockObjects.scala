package cromwell.backend.impl.jes

import cromwell.backend.impl.jes.authentication.JesAuths
import cromwell.filesystems.gcs.auth.GoogleAuthMode

object MockObjects {
  val mockCredentials = JesAuths(GoogleAuthMode.NoAuthMode, GoogleAuthMode.NoAuthMode)
}
