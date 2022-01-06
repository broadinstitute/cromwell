package cromwell.backend.impl.tes

import spray.json._

final case class CreateTaskResponse(id: String)
final case class MinimalTaskView(id: String, state: String)
final case class CancelTaskResponse()

object TesResponseJsonFormatter extends DefaultJsonProtocol {
  implicit val resourcesFormat = jsonFormat6(Resources)
  implicit val inputFormat = jsonFormat6(Input)
  implicit val outputFormat = jsonFormat5(Output)
  implicit val executorFormat = jsonFormat7(Executor)
  implicit val executorLogFormat = jsonFormat5(ExecutorLog)
  implicit val outputFileLogFormat = jsonFormat3(OutputFileLog)
  implicit val taskLogFormat = jsonFormat6(TaskLog)
  implicit val taskFormat = jsonFormat11(Task)
  implicit val minimalTaskView = jsonFormat2(MinimalTaskView)
  implicit val createTaskResponseFormat = jsonFormat1(CreateTaskResponse)
  implicit val cancelTaskResponseFormat = jsonFormat0(CancelTaskResponse)
}
