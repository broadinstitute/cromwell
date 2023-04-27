package cromwell.backend.google.pipelines.batch.authentication

import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

case class GcpBatchAuths(genomics: GoogleAuthMode, gcs: GoogleAuthMode)
