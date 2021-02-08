package cromwell.backend.google.pipelines.common.action

object ActionLabels {
  object Key {
    /**
      * Very short description of the action
      */
    val Tag = "tag"
    val Logging = "logging"
    val InputName = "inputName"
    val OutputName = "outputName"
  }
  object Value {
    val ContainerSetup = "ContainerSetup"
    val UserAction = "UserAction"
    val Localization = "Localization"
    val Delocalization = "Delocalization"
    val Monitoring = "Monitoring"
    val Background = "Background"
    val RetryWithMoreMemory = "CheckingForMemoryRetry"
  }
}
