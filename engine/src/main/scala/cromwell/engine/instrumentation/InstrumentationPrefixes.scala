package cromwell.engine.instrumentation

object InstrumentationPrefixes {
  val WorkflowPrefix: Option[String] = Option("workflow")
  val JobPrefix: Option[String] = Option("job")
  val IoPrefix: Option[String] = Option("io")
  val ApiPrefix: Option[String] = Option("http")
}
