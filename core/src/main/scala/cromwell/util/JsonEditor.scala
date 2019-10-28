package cromwell.util

import cats.data.NonEmptyList
import common.collections.EnhancedCollections._
import cromwell.core.WorkflowId
import io.circe.{Json, JsonNumber, JsonObject, Printer}
import mouse.all._
import io.circe.Json.Folder

import scala.collection.immutable
import common.util.StringUtil._

object JsonEditor {

  private val subWorkflowMetadataKey = "subWorkflowMetadata"
  private val subWorkflowIdKey = "subWorkflowId"

  def includeExcludeJson(json: Json, includeKeys: Option[NonEmptyList[String]], excludeKeys: Option[NonEmptyList[String]]): Json =
    (includeKeys, excludeKeys) match {
      // Take includes, and then remove excludes
      case (Some(includeKeys), Some(excludeKeys)) => includeJson(json, includeKeys) |> (excludeJson(_, excludeKeys))
      case (None, Some(excludeKeys)) => excludeJson(json, excludeKeys)
      case (Some(includeKeys), None) => includeJson(json, includeKeys)
      case _ => json
    }

  def includeJson(json: Json, keys: NonEmptyList[String]): Json = {
    val keysWithId = "id" :: keys
    def folder: Folder[(Json, Boolean)] = new Folder[(Json, Boolean)] {
      override def onNull: (Json, Boolean) = (Json.Null, false)
      override def onBoolean(value: Boolean): (Json, Boolean) = (Json.fromBoolean(value), false)
      override def onNumber(value: JsonNumber): (Json, Boolean) = (Json.fromJsonNumber(value), false)
      override def onString(value: String): (Json, Boolean) = (Json.fromString(value), false)
      override def onArray(value: Vector[Json]): (Json, Boolean) = {
        val newArrayAndKeeps: immutable.Seq[(Json, Boolean)] = value.map(_.foldWith(folder))
        val keep: Boolean = newArrayAndKeeps.map{ case (_, keep) => keep}.foldLeft(false)(_ || _)
        (Json.fromValues(newArrayAndKeeps.map{ case (newJson, _) => newJson}), keep)
      }

      override def onObject(value: JsonObject): (Json, Boolean) = {
        val modified: immutable.List[(String, Json)] = value.toList.flatMap{
          case (key, value) =>
            val keep = keysWithId.foldLeft(false)(_ || key.startsWith(_))
            if (keep)
              List[(String,Json)]((key,value))
            else {
              //run against children, if none of the children need it we can throw it away
              val newJsonAndKeep: (Json, Boolean) = value.foldWith(folder)
              val (newJson, keepChildren) = newJsonAndKeep
              if (keepChildren)
                List((key,newJson))
              else
                List.empty[(String,Json)]
            }
        }
        val jsonObject = Json.fromJsonObject(JsonObject.fromIterable(modified))
        val keep = modified.nonEmpty
        (jsonObject, keep)
      }
    }
    val (newJson,_) = json.foldWith(folder)
    newJson
  }

  def excludeJson(json: Json, keys: NonEmptyList[String]): Json = {
    json.withObject{obj =>
      val keysFiltered = obj.filterKeys(key => !keys.foldLeft(false)(_ || key.startsWith(_)))
      val childrenMapped = keysFiltered.mapValues(excludeJson(_, keys))
      Json.fromJsonObject(childrenMapped)
    }.withArray{
      array =>
        val newArray = array.map(excludeJson(_, keys))
        Json.fromValues(newArray)
    }
  }

  def extractCall(json: Json, callFqn: String, index: Option[Int], attempt: Option[Int]): Json = json // ??? // FIXME

  def removeSubworkflowMetadata(json: Json): Json = json // ??? // FIXME

  def outputs(json: Json): Json = includeJson(json, NonEmptyList.of("outputs")) |> (excludeJson(_, NonEmptyList.one("calls")))

  def logs(json: Json): Json = includeJson(json, NonEmptyList.of("stdout", "stderr", "backendLogs"))

  implicit class EnhancedJson(val json: Json) extends AnyVal {
    def workflowId: Either[Exception, WorkflowId] = {
      val workflowIdOpt = for {
        o <- json.asObject
        id <- o.kleisli("id")
        s <- id.asString
      } yield WorkflowId.fromString(s)

      workflowIdOpt match {
        case Some(workflowId) => Right(workflowId)
        case None => Left(new RuntimeException(s"did not find workflow id in ${json.printWith(Printer.spaces2).elided(100)}"))
      }
    }
  }

  /**
    * Look for an optional JsonObject by its key
    * @param workflowJson - workflow Json to look in
    * @param key - key to look for
    * @return - optional tuple of workflow JsonObject and found element JsonObject
    */
  private def extractJsonObjectByKey(workflowJson: Json, key: String): Option[(JsonObject, JsonObject)] =
    extractJsonByKey(workflowJson, key) flatMap {
      case (jsonObj, json) => json.asObject.flatMap(jo => Option((jsonObj, jo)))
    }

  /**
    * Look for an optional JsonObject by its key
    * @param workflowJson - workflow Json to look in
    * @param key - key to look for
    * @return - optional tuple of workflow JsonObject and found element Json
    */
  private def extractJsonByKey(workflowJson: Json, key: String): Option[(JsonObject, Json)] =
    for {
      wo <- workflowJson.asObject
      foundJson <- wo(key)
    } yield (wo, foundJson)

  /**
    * In-memory upsert of labels into the base Json, handling root and sub workflows appropriately.
    *
    * @param json json blob with or without "labels" field
    * @param databaseLabels a map of workflow IDs to maps of labels one would like to apply to a workflow json
    * @return json with labels merged in.  Any prior non-object "labels" field will be overwritten and any object fields will be merged together and - again - any existing values overwritten.
    */
  def updateLabels(json: Json, databaseLabels: Map[WorkflowId, Map[String, String]]): Json = {

    def updateLabelsInCalls(callObject: JsonObject, subworkflowJson: Json): Json = {
      // If the call contains a subWorkflowMetadata key, return a copy of the call with
      // its subworkflowMetadata updated via a recursive call to `doUpdateWorkflow`.
      val updatedSubworkflow = doUpdateWorkflow(subworkflowJson)
      Json.fromJsonObject(callObject.add(subWorkflowMetadataKey, updatedSubworkflow))
    }

    def doUpdateWorkflow(workflowJson: Json): Json = {
      val id: String = workflowJson.workflowId match {
        case Right(id) => id.toString
        case Left(e) => throw e
      }

      val workflowWithUpdatedCalls = updateWorkflowCallsJson(workflowJson, updateLabelsInCalls)

      databaseLabels.get(WorkflowId.fromString(id)) match {
        case None => workflowWithUpdatedCalls
        case Some(labels) =>
          val labelsJson: Json = Json.fromFields(labels safeMapValues Json.fromString)
          workflowWithUpdatedCalls deepMerge Json.fromFields(List(("labels", labelsJson)))
      }
    }

    doUpdateWorkflow(workflowJson = json)
  }

  /**
    * Update workflow json, replacing "subWorkflowMetadata" elements of root workflow's "calls" object by "subworkflowId"
    *
    * @param workflowJson json blob
    * @return updated json
    */
  def replaceSubworkflowMetadataWithId(workflowJson: Json): Json = {
    def replaceSubworkflowMetadataObjectWithSubworkflowIdInCalls(callObject: JsonObject, subworkflowJson: Json): Json = {
      val subWorkflowId = subworkflowJson.workflowId match {
        case Right(id) => id.toString
        case Left(e) => throw e
      }
      val updatedCallObj = callObject.remove(subWorkflowMetadataKey).add(subWorkflowIdKey, Json.fromString(subWorkflowId))
      Json.fromJsonObject(updatedCallObj)
    }

    updateWorkflowCallsJson(workflowJson, replaceSubworkflowMetadataObjectWithSubworkflowIdInCalls)
  }

  private def updateWorkflowCallsJson(workflowJson: Json, updateCallsFunc: (JsonObject, Json) => Json): Json = {
    val workflowWithUpdatedCalls: Json = extractJsonObjectByKey(workflowJson, "calls") match {
      // If there were no calls just return the workflow JSON unmodified.
      case None => workflowJson
      case Some((_, calls)) =>
        val updatedCalls = calls.mapValues {
          // The Json (a JSON array, really) corresponding to the array of call objects for a call name.
          callValue: Json =>
            // The object above converted to a List[Json].
            val callArray: Vector[Json] = callValue.asArray.toVector.flatten

            val updatedCallArray = callArray map { callJson =>
              // If there is no subworkflow object this will be None.
              val callAndSubworkflowObjects: Option[(JsonObject, Json)] = extractJsonByKey(callJson, subWorkflowMetadataKey)

              callAndSubworkflowObjects match {
                case None => callJson
                case Some((callObject, subworkflowJson)) => updateCallsFunc(callObject, subworkflowJson)
              }
            }
            Json.fromValues(updatedCallArray)
        }
        Json.fromJsonObject(workflowJson.asObject.get.add("calls", Json.fromJsonObject(updatedCalls)))
    }
    workflowWithUpdatedCalls
  }
}
