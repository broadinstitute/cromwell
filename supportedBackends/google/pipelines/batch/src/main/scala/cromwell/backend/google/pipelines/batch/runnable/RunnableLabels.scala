package cromwell.backend.google.pipelines.batch.runnable

object RunnableLabels {
  object Key {
    /**
      * Very short description of the runnable
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
