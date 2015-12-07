package cromwell.webservice

import cromwell.binding.FullyQualifiedName
import cromwell.binding.values.WdlValue
import cromwell.engine.backend.{WorkflowQueryResult, CallMetadata, CallLogs}
import org.joda.time.DateTime
import spray.json.JsObject

case class WorkflowValidateResponse(valid: Boolean, error: Option[String])

case class WorkflowStatusResponse(id: String, status: String)

case class WorkflowSubmitResponse(id: String, status: String)

case class WorkflowOutputResponse(id: String, outputs: Map[FullyQualifiedName, WdlValue])

case class WorkflowAbortResponse(id: String, status: String)

case class CallOutputResponse(id: String, callFqn: String, outputs: Map[FullyQualifiedName, WdlValue])

case class CallStdoutStderrResponse(id: String, logs: Map[String, Seq[CallLogs]])

case class WorkflowMetadataResponse(id: String, status: String, submission: DateTime, start: Option[DateTime],
                                    end: Option[DateTime], inputs: JsObject, outputs: Option[Map[String, WdlValue]],
                                    calls: Map[String, Seq[CallMetadata]])

case class WorkflowQueryResponse(results: Seq[WorkflowQueryResult])

final case class CallCachingResponse(updateCount: Int)
