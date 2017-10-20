package cromwell.backend.impl.tes

import spray.json._

final case class CreateTaskResponse(id: String)
final case class MinimalTaskView(id: String, state: String)
final case class CancelTaskResponse()

object TesResponseJsonFormatter extends DefaultJsonProtocol {
  implicit val resourcesFormat = jsonFormat5(Resources)
  implicit val taskParameterFormat = jsonFormat6(TaskParameter)
  implicit val portsFormat = jsonFormat2(Ports)
  implicit val executorFormat = jsonFormat8(Executor)
  implicit val executorLogFormat = jsonFormat7(ExecutorLog)
  implicit val outputFileLogFormat = jsonFormat3(OutputFileLog)
  implicit val taskLogFormat = jsonFormat5(TaskLog)
  implicit val taskFormat = jsonFormat12(Task)
  implicit val minimalTaskView = jsonFormat2(MinimalTaskView)
  implicit val createTaskResponseFormat = jsonFormat1(CreateTaskResponse)
  implicit val cancelTaskResponseFormat = jsonFormat0(CancelTaskResponse)
}
