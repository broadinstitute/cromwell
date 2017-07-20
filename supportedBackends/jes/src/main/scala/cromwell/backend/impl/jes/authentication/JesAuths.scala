package cromwell.backend.impl.jes.authentication

import cromwell.cloudSupport.gcp.auth.GoogleAuthMode

case class JesAuths(genomics: GoogleAuthMode, gcs: GoogleAuthMode)
