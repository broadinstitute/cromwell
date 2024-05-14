package cromwell.services.metrics.bard.model

case class TaskSummaryEvent(terminationStatus: String, returnCode: Int) extends BardEvent {
  override def eventName: String = "task:summary"

}
