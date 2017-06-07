package cromwell.webservice

import java.nio.file.Paths
import java.time.OffsetDateTime

import cromwell.core._
import cromwell.engine._
import cromwell.services.metadata.MetadataService
import MetadataService.{WorkflowQueryResponse, WorkflowQueryResult}
import cromwell.util.JsonFormatting.WdlValueJsonFormatter
import WdlValueJsonFormatter._
import better.files.File
import cromwell.database.sql.tables.CallCachingEntry
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor.{CallCachingDiff, CallCachingDiffElement}
import spray.json._

object WorkflowJsonSupport extends DefaultJsonProtocol {
  implicit val workflowStatusResponseProtocol = jsonFormat2(WorkflowStatusResponse)
  implicit val workflowAbortResponseProtocol = jsonFormat2(WorkflowAbortResponse)
  implicit val workflowSubmitResponseProtocol = jsonFormat2(WorkflowSubmitResponse)
  implicit val workflowOutputResponseProtocol = jsonFormat2(WorkflowOutputResponse)
  implicit val callOutputResponseProtocol = jsonFormat3(CallOutputResponse)
  implicit val engineStatsProtocol = jsonFormat2(EngineStatsActor.EngineStats)
  implicit val callAttempt = jsonFormat2(CallAttempt)
  implicit val workflowSourceData = jsonFormat4(WorkflowSourceFilesWithoutImports)

  implicit object fileJsonFormat extends RootJsonFormat[File] {
    override def write(obj: File) = JsString(obj.path.toAbsolutePath.toString)
    override def read(json: JsValue): File = json match {
      case JsString(str) => Paths.get(str)
      case unknown => throw new NotImplementedError(s"Cannot parse $unknown to a File")
    }
  }


  implicit object callCachingDiffElementJsonFormat extends RootJsonFormat[CallCachingDiffElement] {
    override def write(obj: CallCachingDiffElement) =
      JsObject(Map(
        obj.hashKey -> JsObject(Map(
          "callA" -> obj.hashValueA.toJson,
          "callB" -> obj.hashValueB.toJson
        ))
      ))
    override def read(json: JsValue): CallCachingDiffElement =
      throw new NotImplementedError(s"Cannot parse json to CallCachingDiffElement")
  }
  
  implicit object callCachingDiffJsonFormat extends RootJsonFormat[CallCachingDiff] {
    def makeCallObject(cacheEntry: CallCachingEntry) =  JsObject(Map(
      "workflowId" -> JsString(cacheEntry.workflowExecutionUuid),
      "callFqn" -> JsString(cacheEntry.callFullyQualifiedName),
      "jobIndex" -> JsNumber(cacheEntry.jobIndex),
      "allowResultReuse" -> JsBoolean(cacheEntry.allowResultReuse)
    ))
    
    override def write(obj: CallCachingDiff) = 
      JsObject(Map(
        "callA" -> makeCallObject(obj.cacheEntryA),
        "callB" -> makeCallObject(obj.cacheEntryB),
        "hashDifferential" -> obj.hashDifferential.toList.toJson
      ))
    override def read(json: JsValue): CallCachingDiff = 
      throw new NotImplementedError(s"Cannot parse json to CallCachingDiff")
  }

  implicit val workflowSourceDataWithImports = jsonFormat5(WorkflowSourceFilesWithDependenciesZip)
  implicit val errorResponse = jsonFormat3(FailureResponse)
  implicit val successResponse = jsonFormat3(SuccessResponse)

  implicit object DateJsonFormat extends RootJsonFormat[OffsetDateTime] {
    override def write(obj: OffsetDateTime) = JsString(obj.toString)

    override def read(json: JsValue): OffsetDateTime = json match {
      case JsString(str) => OffsetDateTime.parse(str)
      case unknown => throw new NotImplementedError(s"Cannot parse $unknown to a DateTime")
    }
  }

  implicit val workflowQueryResult = jsonFormat5(WorkflowQueryResult)
  implicit val workflowQueryResponse = jsonFormat1(WorkflowQueryResponse)
}

