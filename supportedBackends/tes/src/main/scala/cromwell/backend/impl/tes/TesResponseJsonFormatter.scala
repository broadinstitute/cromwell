package cromwell.backend.impl.tes

import spray.json._

final case class CreateTaskResponse(id: String)
final case class MinimalTaskView(id: String, state: String)
final case class CancelTaskResponse()

object TesResponseJsonFormatter extends DefaultJsonProtocol {

  /**Spray-json serializes Scala `None` as JSON `null`.
   * Protobuf3's JSON parser (eg in Funnel/TES) rejects null for string/message fields:
   *   "invalid value for string field value: null"
   * Omitting the key entirely is the correct proto3-compliant behaviour — missing
   * fields are treated as their default value (empty string / zero / etc.).
   *
   * Two sources of null in the TES payload:
   *   1. Top-level Option fields on Task/Input/Output/etc. (e.g. Task.id, Input.content)
   *   2. Map[String, Option[String]] values inside `tags` and `backend_parameters`
   *      (e.g. Task.tags["parent_workflow_id"] = None for top-level workflows)
   *
   * We use a deep recursive strip so that nulls are removed at every nesting level.
   * Reading is unaffected: spray-json already maps missing JSON keys to `None`.
   */
  private def deepStripNulls(json: JsValue): JsValue = json match {
    case o: JsObject =>
      JsObject(o.fields.filterNot(_._2 == JsNull).map { case (k, v) => k -> deepStripNulls(v) })
    case a: JsArray =>
      JsArray(a.elements.filterNot(_ == JsNull).map(deepStripNulls))
    case other => other
  }

  private def nullFreeFormat[T](fmt: RootJsonFormat[T]): RootJsonFormat[T] =
    new RootJsonFormat[T] {
      def read(json: JsValue): T = fmt.read(json)
      def write(obj: T): JsValue = deepStripNulls(fmt.write(obj))
    }

  /** OutputFileLog is partly defined by size_bytes: Long which is supplied as string as described by the TES spec.
   * We create a custom deserializer to read the size_bytes string and turn it into the Long Cromwell expects.
   * NOTE: size_bytes must be Long, not Int — output files can exceed 2 GB (Int.MaxValue = 2,147,483,647).
   */
  implicit object customJsonFormatOutputFileLog extends RootJsonFormat[OutputFileLog] {
    def write(obj: OutputFileLog): JsValue =
      JsObject(
        "url" -> JsString(obj.url),
        "path" -> JsString(obj.path),
        "size_bytes" -> JsNumber(obj.size_bytes)
      )

    def read(value: JsValue): OutputFileLog =
      value.asJsObject.getFields("url", "path", "size_bytes") match {
        case Seq(JsString(url), JsString(path), JsString(size_bytes)) =>
          OutputFileLog(url, path, size_bytes.toLong)
        case Seq(JsString(url), JsString(path)) => OutputFileLog(url, path, 0L)
        case _ => throw DeserializationException("Cannot deserialize OutputFileLog")
      }
  }

  implicit val resourcesFormat          = nullFreeFormat(jsonFormat6(Resources))
  implicit val inputFormat              = nullFreeFormat(jsonFormat6(Input))
  implicit val outputFormat             = nullFreeFormat(jsonFormat5(Output))
  implicit val executorFormat           = nullFreeFormat(jsonFormat7(Executor))
  implicit val executorLogFormat        = nullFreeFormat(jsonFormat5(ExecutorLog))
  implicit val taskLogFormat            = nullFreeFormat(jsonFormat6(TaskLog))
  implicit val taskFormat               = nullFreeFormat(jsonFormat11(Task))
  implicit val minimalTaskView          = jsonFormat2(MinimalTaskView)
  implicit val createTaskResponseFormat = jsonFormat1(CreateTaskResponse)
  implicit val cancelTaskResponseFormat = jsonFormat0(CancelTaskResponse)
}
