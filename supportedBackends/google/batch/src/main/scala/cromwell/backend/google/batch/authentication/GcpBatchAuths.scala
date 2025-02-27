package cromwell.backend.google.batch.authentication

import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

case class GcpBatchAuths(batch: GoogleAuthMode, gcs: GoogleAuthMode)
