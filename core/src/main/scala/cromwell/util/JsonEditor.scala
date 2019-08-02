package cromwell.util

import cats.data.NonEmptyList
import common.collections.EnhancedCollections._
import io.circe.{Json, JsonNumber, JsonObject}
import mouse.all._
import io.circe.Json.Folder

import scala.collection.immutable

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
        val keep = modified.size > 0
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

  /**
    * We only return labels associated with top-level workflows.  Subworkflows don't include labels (as of 7/26/19).
    *
    * Thus this method puts the labels as a top-level field.
    *
    * @param json json blob with or without "labels" field
    * @param labels a map of labels one would like to apply to a workflow json
    * @return json with labels merged in.  Any prior non-object "labels" field will be overwritten and any object fields will be merged together and - again - any existing values overwritten.
    */
  def augmentLabels(json: Json, labels: Map[String, String]): Json = {
    val newData: Json = Json.fromFields(labels.safeMapValues(Json.fromString))
    val newObj: Json = Json.fromFields(List(("labels", newData)))
    //in the event of a key clash, the values in "newObj" will be favored over "json"
    json deepMerge newObj
  }

  def removeSubworkflowData(json: Json): Json = excludeJson(json, NonEmptyList.of("subWorkflowMetadata"))
}
