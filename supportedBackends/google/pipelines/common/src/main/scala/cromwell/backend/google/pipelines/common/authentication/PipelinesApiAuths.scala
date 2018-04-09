package cromwell.backend.google.pipelines.common.authentication

import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

case class PipelinesApiAuths(genomics: GoogleAuthMode, gcs: GoogleAuthMode)
