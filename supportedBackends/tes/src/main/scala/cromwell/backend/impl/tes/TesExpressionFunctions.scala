package cromwell.backend.impl.tes

import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.standard.StandardExpressionFunctionsParams

class TesExpressionFunctions(standardParams: StandardExpressionFunctionsParams)
    extends SharedFileSystemExpressionFunctions(standardParams) {

  // The original check only passed through local absolute paths ("/...")
  // and FTP URIs, sending everything else through callContext.root.resolve(str). For cloud-backed TES
  // deployments (S3, GCS, Azure, …) the output paths returned by the TES server are already fully-
  // qualified URIs (e.g. "s3://bucket/key"). Resolving a URI-scheme string against an S3Path doubles
  // the path (s3://…execution/s3:/…execution/stdout) because S3Path.resolve(String) treats any string
  // that does not start with "/" as relative.
  // Fix: treat any string that contains "://" as already-absolute and return it unchanged.
  override def preMapping(str: String) =
    if (str.startsWith("/") || str.contains("://")) str
    else standardParams.callContext.root.resolve(str).toString
}
