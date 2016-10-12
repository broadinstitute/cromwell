package cromwell.backend.impl.jes.authentication

import cromwell.filesystems.gcs.auth.GoogleAuthMode

case class JesAuths(genomics: GoogleAuthMode, gcs: GoogleAuthMode)
