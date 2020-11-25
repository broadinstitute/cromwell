package cromwell.backend.google.pipelines.v2alpha1.api

//noinspection TypeAnnotation
object ActionFlag extends Enumeration {
  type ActionFlag = Value

  val FlagUnspecified = Value("FLAG_UNSPECIFIED")
  val IgnoreExitStatus = Value("IGNORE_EXIT_STATUS")
  val RunInBackground = Value("RUN_IN_BACKGROUND")
  val AlwaysRun = Value("ALWAYS_RUN")
  val EnableFuse = Value("ENABLE_FUSE")
  val PublishExposedPorts = Value("PUBLISH_EXPOSED_PORTS")
  val DisableImagePrefetch = Value("DISABLE_IMAGE_PREFETCH")
}
