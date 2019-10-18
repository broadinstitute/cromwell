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

  def includeExcludeJson(json: Json, includeKeys: Option[NonEmptyList[String]], excludeKeys: Option[NonEmptyList[String]]): Json =
    (includeKeys, excludeKeys) match {
      // Take includes, and then remove excludes
      case (Some(includeKeys), Some(excludeKeys)) => includeJson(json, includeKeys) |> (excludeJson(_, excludeKeys))
      case (None, Some(excludeKeys)) => excludeJson(json, excludeKeys)
      case (Some(includeKeys), None) => includeJson(json, includeKeys)
      case _ => json
    }

  def includeJson(json: Json, keys: NonEmptyList[String]): Json = {
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
            val keep =  keys.foldLeft(false)(_ || key.startsWith(_))
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

  def outputs(json: Json): Json = includeJson(json, NonEmptyList.of("outputs", "id")) |> (excludeJson(_, NonEmptyList.one("calls")))

  def logs(json: Json): Json = includeJson(json, NonEmptyList.of("stdout", "stderr", "backendLogs", "id"))

  implicit class EnhancedJson(val json: Json) extends AnyVal {
    def rootWorkflowId: Option[WorkflowId] = {
      for {
        o <- json.asObject
        id <- o.kleisli("id")
        s <- id.asString
      } yield WorkflowId.fromString(s)
    }
  }

  /**
    * In-memory upsert of labels into the base Json, handling root and sub workflows appropriately.
    *
    * @param json json blob with or without "labels" field
    * @param databaseLabels a map of workflow IDs to maps of labels one would like to apply to a workflow json
    * @return json with labels merged in.  Any prior non-object "labels" field will be overwritten and any object fields will be merged together and - again - any existing values overwritten.
    */
  def updateLabels(json: Json, databaseLabels: Map[WorkflowId, Map[String, String]]): Json = {
    val subWorkflowMetadataKey = "subWorkflowMetadata"

    def doUpdateWorkflow(workflowJson: Json): Json = {
      val id: String = workflowJson.rootWorkflowId.
        getOrElse(throw new RuntimeException(s"did not find workflow id in ${workflowJson.printWith(Printer.spaces2).elided(100)}")).toString

      // Look for an optional JsonObject keyed by "calls".
      val callsObject: Option[JsonObject] = for {
        wo <- workflowJson.asObject
        callsJson <- wo("calls")
        co <- callsJson.asObject
      } yield co

      val workflowWithUpdatedCalls: Json = callsObject match {
        // If there were no calls just return the workflow JSON unmodified.
        case None => workflowJson
        case Some(calls) =>
          val updatedCallsObject = calls.mapValues {
            // The Json (a JSON array, really) corresponding to the array of call objects for a call name.
            callValue: Json =>
              // The object above converted to a List[Json].
              val callArray: Vector[Json] = callValue.asArray.toVector.flatten

              val updatedCallArray = callArray map { callJson =>
                // If there is no subworkflow object this will be None.
                val callAndSubworkflowObjects: Option[(JsonObject, Json)] = for {
                  co <- callJson.asObject
                  sub <- co(subWorkflowMetadataKey)
                  _ <- sub.asObject
                } yield (co, sub)

                callAndSubworkflowObjects match {
                  case None => callJson
                  case Some((callObject, subworkflowObject)) =>
                    // If the call contains a subWorkflowMetadata key, return a copy of the call with
                    // its subworkflowMetadata updated via a recursive call to `doUpdateWorkflow`.
                    val updatedSubworkflow = doUpdateWorkflow(subworkflowObject)
                    Json.fromJsonObject(callObject.add(subWorkflowMetadataKey, updatedSubworkflow))
                }
              }
              Json.fromValues(updatedCallArray)
          }
          Json.fromJsonObject(workflowJson.asObject.get.add("calls", Json.fromJsonObject(updatedCallsObject)))
      }

      databaseLabels.get(WorkflowId.fromString(id)) match {
        case None => workflowWithUpdatedCalls
        case Some(labels) =>
          val labelsJson: Json = Json.fromFields(labels safeMapValues Json.fromString)
          workflowWithUpdatedCalls deepMerge Json.fromFields(List(("labels", labelsJson)))
      }
    }

    doUpdateWorkflow(workflowJson = json)
  }
}
