package cromwell.backend.google.batch.models

case class ManifestFile(imageIdentifier: String, diskSizeGb: Int, files: List[ReferenceFile])
