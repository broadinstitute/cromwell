package cromwell.core.instrumentation

object InstrumentationPrefixes {
  val ApiPrefix: Option[String] = Option("rest-api")
  val BackendPrefix: Option[String] = Option("backend")
  val JobPrefix: Option[String] = Option("job")
  val IoPrefix: Option[String] = Option("io")
  val ServicesPrefix: Option[String] = Option("services")
  val WorkflowPrefix: Option[String] = Option("workflow")
  val CromIamSamOverheadPrefix: Option[String] = Option("sam-overhead")
}
