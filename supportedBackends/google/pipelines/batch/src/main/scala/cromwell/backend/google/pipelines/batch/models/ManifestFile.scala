package cromwell.backend.google.pipelines.batch.models

case class ManifestFile(imageIdentifier: String, diskSizeGb: Int, files: List[ReferenceFile])
