package cromwell.backend.impl.jes.authentication

import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

case class JesAuths(genomics: GoogleAuthMode, gcs: GoogleAuthMode)
