package cromwell.backend.impl.tes

import spray.json.{RootJsonFormat, _}

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
//  implicit val taskLogFormat = jsonFormat6(TaskLog)
//  implicit val taskFormat = jsonFormat11(Task)
//  implicit val minimalTaskView = jsonFormat2(MinimalTaskView)
  implicit val createTaskResponseFormat = jsonFormat1(CreateTaskResponse)
  implicit val cancelTaskResponseFormat = jsonFormat0(CancelTaskResponse)

  implicit object myCustomJsonFormatMinimal extends RootJsonFormat[MinimalTaskView] {
    def read(value: JsValue) = {
      System.out.print("CUSTOM FORMATTER MINIMAL    " + value)
      jsonFormat2(MinimalTaskView).read(value)
    }
    def write(obj: MinimalTaskView): JsValue = jsonFormat2(MinimalTaskView).write(obj)
  }

  implicit object myCustomJsonFormatTaskLog extends RootJsonFormat[TaskLog] {
    def read(value: JsValue) = {
      System.out.print("CUSTOM FORMATTER TASKLOG      " + value)
      jsonFormat6(TaskLog).read(value)
    }
    def write(obj: TaskLog): JsValue = jsonFormat6(TaskLog).write(obj)
  }

  implicit object myCustomJsonFormatTask extends RootJsonFormat[Task] {
    def read(value: JsValue) = {
      System.out.print("CUSTOM FORMATTER TASK     " + value)
      jsonFormat11(Task).read(value)
    }
    def write(obj: Task): JsValue = jsonFormat11(Task).write(obj)
  }
}
