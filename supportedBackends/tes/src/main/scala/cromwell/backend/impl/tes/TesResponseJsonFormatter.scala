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
  implicit val taskLogFormat = jsonFormat6(TaskLog)
  implicit val taskFormat = jsonFormat11(Task)
  implicit val minimalTaskView = jsonFormat2(MinimalTaskView)
  implicit val createTaskResponseFormat = jsonFormat1(CreateTaskResponse)
  implicit val cancelTaskResponseFormat = jsonFormat0(CancelTaskResponse)

  /** OutputFileLog is partly defined by size_byte: Int which is supplied as string as described by the TES spec.
   * We create a custom deserializer to read the size_byte string and turn it into the Int Cromwell expects
   */
  implicit object customJsonFormatOutputFileLog extends RootJsonFormat[OutputFileLog] {
    def write(obj: OutputFileLog): JsValue = JsArray(JsString(obj.url), JsString(obj.path), JsNumber(obj.size_bytes))

    def read(value: JsValue) = {
      System.out.print("READ OUTPUTFILELOG:       " + value)
      jsonFormat3(OutputFileLog).read(value)
    } // value match {
//      case JsArray(Vector(JsString(url), JsString(path), JsString(size_bytes))) =>
//        OutputFileLog(url, path, size_bytes.toInt)
//      case _ => deserializationError("Cannot deserialize into OutputFileLog")
//    }
  }
}
