package cromwell.util

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.instances.vector._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.collections.EnhancedCollections._
import common.util.StringUtil._
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.core.WorkflowId
import io.circe.Json.Folder
import io.circe.{Json, JsonNumber, JsonObject, Printer}

import scala.collection.immutable

object JsonEditor {

  private val subWorkflowMetadataKey = "subWorkflowMetadata"
  private val subWorkflowIdKey = "subWorkflowId"
  private val keysToIncludeInCallsOrWorkflows = NonEmptyList.of("id", "shardIndex", "attempt")

  def includeExcludeJson(json: Json, includeKeys: Option[NonEmptyList[String]], excludeKeys: Option[NonEmptyList[String]]): ErrorOr[Json] =
    (includeKeys, excludeKeys) match {
      // Take includes, and then remove excludes
      case (Some(includeKeys), Some(excludeKeys)) => includeJson(json, includeKeys) flatMap (excludeJson(_, excludeKeys))
      case (None, Some(excludeKeys)) => excludeJson(json, excludeKeys)
      case (Some(includeKeys), None) => includeJson(json, includeKeys)
      case _ => json.validNel
    }

  def includeJson(json: Json, keys: NonEmptyList[String]): ErrorOr[Json] = {
    val keysWithId = keysToIncludeInCallsOrWorkflows ::: keys
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
            val keep = keysWithId.foldLeft(false)(_ || key.equals(_))
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
    newJson.validNel
  }

  case class Filter(components: NonEmptyList[String]) {
    def isSameOrMoreSpecificThan(that: Filter): Boolean = {
      // A request to remove a filter of 'outputs' should also remove a more specific 'outputs:foo' filter if present.
      // This algorithm effectively converts `Filter`s to representations in classic metadata syntax to see if `that`
      // filter is a prefix of `this` filter, i.e. if `this` is a more specific filter of `that`.
      def toClassic(filter: Filter): String = filter.components.toList.mkString("", ":", ":")
      val List(thisComponents, thatComponents) = List(this, that) map toClassic
      thisComponents startsWith thatComponents
    }
  }

  object Filter {
    def fromString(string: String): Filter = {
      string.split(':').toList match {
        case Nil => throw new RuntimeException("Programmer error: unpossible Nil result from split")
        case h :: t => Filter(NonEmptyList.of(h, t: _*))
      }
    }
  }

  sealed trait FilterResult
  case object Keep extends FilterResult
  case object Discard extends FilterResult
  case class Replace(updatedJson: Json) extends FilterResult

  /** A `FilterGroup` represents all the `Filter`s for include xor exclude. The argument is intentionally not a NEL
    * since `FilterGroup`s should start out non-empty but then become empty after `remove`s. e.g. excludeKey 'id'. */
  case class FilterGroup(filters: List[Filter]) {

    def remove(head: String, tail: String*): FilterGroup = {
      remove(Filter.fromString(head), tail map Filter.fromString: _*)
    }

    /** Return a copy of this `FilterGroup` with all of the specified filters removed. */
    def remove(head: Filter, tail: Filter*): FilterGroup = {
      (head :: tail.toList).foldLeft(this) { case (acc, f) =>
        // Copy the `acc` `FilterGroup` with any filters that are the same or more specific than `f` removed.
        acc.copy(acc.filters.filterNot { _.isSameOrMoreSpecificThan(f) })
      }
    }

    def filter(jsonObject: JsonObject): JsonObject = {
      // If at the top level there is a perfect match, discard.
      // If at the top level there are no matches at all, keep.
      // Otherwise fold.
      def filtersMatchingKey(objectKey: String): List[Filter] = filters.filter(_.components.head == objectKey)

      def descend(json: Json, remainingComponents: NonEmptyList[String]): Json = {
        json.asObject match {
          case None => json
          case Some(jsonObject) =>
            val component = remainingComponents.head
            jsonObject(component) match {
              case None => json
              case Some(inner) =>
                val updatedObject = remainingComponents.tail match {
                  case Nil => jsonObject.remove(component)
                  case h :: t => jsonObject.add(component, descend(inner, NonEmptyList.of(h, t: _*)))
                }
                Json.fromJsonObject(updatedObject)
            }
        }
      }

      val updatedKeyValues: List[(String, Json)] = jsonObject.toList flatMap { case (key, json) =>
        filtersMatchingKey(key) match {
          case Nil =>
            // No filters matched this object key, return the key/value unmodified.
            List((key, json))
          case fs if fs.exists(_.components.tail == Nil) =>
            // If there is a filter that has no other components and thus matches the object key completely,
            // return a Nil List to delete the key/value pair from the object.
            Nil
          case fs =>
            // The JSON value will not be deleted but may need to be edited. Fold the JSON through all the
            // `fs` filters.
            val possiblyUpdatedJson = fs.foldLeft(json) { case (j, f) =>
              f.components.tail match {
                // We know all of these filters did not have have "tails", i.e. did not
                // match the object key and thus signal that the key/value should be removed.
                case Nil => throw new RuntimeException("Programmer error: filter components tail should not be empty")
                case c :: cs => descend(j, NonEmptyList.of(c, cs: _*))
              }
            }
            List((key, possiblyUpdatedJson))
        }
      }
      JsonObject.fromIterable(updatedKeyValues)
    }
  }

  def excludeJson(json: Json, keys: NonEmptyList[String]): ErrorOr[Json] = {
    val filters: ErrorOr[List[Filter]] = keys.toList traverse { key => key.split(':').toList match {
      case Nil => s"Programmer error: string split resulting in empty array: $key".invalidNel
      case h :: t => Filter(NonEmptyList.of(h, t: _*)).validNel
    }}

    filters map FilterGroup.apply flatMap excludeWorkflowJson(json)
  }

  def excludeWorkflowJson(workflowJson: Json)(filters: FilterGroup): ErrorOr[Json] = {
    // Always include 'id' in workflows, 'shardIndex' and 'attempt' in jobs. i.e. make sure these keys are excluded
    // from the excludes.
    val workflowKeysToExclude = filters.remove("id")
    val jobKeysToExclude = filters.remove("shardIndex", "attempt")

    def excludeJob(job: JsonObject): ErrorOr[JsonObject] = jobKeysToExclude.filter(job).validNel

    def excludeCall(callJson: Json): ErrorOr[Json] = {
      def excludeSubworkflowCall(parentWorkflowCall: JsonObject)(subworkflow: Json): ErrorOr[JsonObject] = {
        excludeWorkflowJson(subworkflow)(filters) map { excludedSubworkflow => parentWorkflowCall.add(subWorkflowMetadataKey, excludedSubworkflow) }
      }

      for {
        callObject <- callJson.asObject.map(_.validNel).getOrElse(s"call JSON unexpectedly not an object: $callJson".invalidNel)
        excludedJob <- excludeJob(callObject)
        excludedCall <- excludedJob(subWorkflowMetadataKey) map excludeSubworkflowCall(excludedJob) getOrElse excludedJob.validNel
      } yield Json.fromJsonObject(excludedCall)
    }

    def excludeJsonArray(callsArrayJson: Json): ErrorOr[Json] = {
      for {
        array <- callsArrayJson.asArray.map(_.validNel).getOrElse(s"calls JSON unexpectedly not an array: $callsArrayJson".invalidNel)
        updatedCalls <- array.traverse[ErrorOr, Json](excludeCall)
      } yield Json.fromValues(updatedCalls)
    }

    def excludeCallsObject(jsonObject: JsonObject): ErrorOr[Json] =
      (jsonObject traverse excludeJsonArray) map Json.fromJsonObject

    for {
      workflowObject <- workflowJson.asObject.map(_.validNel).getOrElse(s"Workflow JSON unexpectedly not an object: $workflowJson".invalidNel)
      workflowObjectWithKeysFiltered = workflowKeysToExclude.filter(workflowObject)
      updatedWorkflowJsonObject <- workflowObjectWithKeysFiltered("calls") match {
        case None => workflowObjectWithKeysFiltered.validNel
        case Some(calls) => calls.asObject match {
          case None => s"calls JSON unexpectedly not an object: $calls".invalidNel
          case Some(callsObject) => excludeCallsObject(callsObject) map {
            updatedCallsObject => workflowObjectWithKeysFiltered.add("calls", updatedCallsObject)
          }
        }
      }
    } yield Json.fromJsonObject(updatedWorkflowJsonObject)
  }

  def outputs(json: Json): ErrorOr[Json] = excludeJson(json, NonEmptyList.one("calls")) flatMap (includeJson(_, NonEmptyList.of("outputs")))

  // Return the calls object with subworkflows removed.
  private def removeSubworkflowCalls(callsObject: JsonObject): ErrorOr[Json] = {
    def isCallSubworkflow(callArrayJson: Json): ErrorOr[Boolean] = {
      for {
        array <- callArrayJson.asArray map (_.validNel) getOrElse s"call array unexpectedly not an array: $callArrayJson".invalidNel
        // All calls within this array are assumed to have the same "shape": either they are subworkflows or regular jobs.
        // So this only looks at the first element of the call array to determine whether this member of the containing
        // object should be filtered as a subworkflow. There should always be at least one element in this calls array.
        callJson <- array.headOption map (_.validNel) getOrElse "call array unexpectedly empty".invalidNel
        call <- callJson.asObject map (_.validNel) getOrElse s"Call JSON unexpectedly not an object: $callJson".invalidNel
      } yield call.contains(subWorkflowMetadataKey) || call.contains(subWorkflowIdKey)
    }

    callsObject.toList.flatTraverse[ErrorOr, (String, Json)] {
      case (key, json) => isCallSubworkflow(json) map { isSubworkflow =>
        // If the value is true return an empty List, otherwise return a single-element List containing the (key, json) pair.
        if (isSubworkflow) List.empty else List((key, json))
      }
    } map { Json.fromFields }
  }

  def logs(workflowJson: Json): ErrorOr[Json] = {
    val inputsAndOutputs = NonEmptyList.of("outputs", "inputs")
    val shardAttemptAndLogsFields = NonEmptyList.of("shardIndex", "attempt", "stdout", "stderr", "backendLogs")

    def updatedWorkflowJson(calls: JsonObject): ErrorOr[Json] = {
      for {
        updatedCallsJson <- removeSubworkflowCalls(calls)
      } yield Json.fromJsonObject(workflowJson.asObject.get.add("calls", updatedCallsJson))
    }

    for {
      calls <- callsObject(workflowJson)
      workflowWithSubworkflowsRemoved <- calls match {
        case None => workflowJson.validNel
        case Some(cs) => updatedWorkflowJson(cs)
      }
      // exclude outputs and inputs since variables can be named anything including internally reserved words like
      // `stdout` and `stderr` which would be erroneously included among the logs.
      excluded <- excludeJson(workflowWithSubworkflowsRemoved, inputsAndOutputs)
      inc <- includeJson(excluded, shardAttemptAndLogsFields)
    } yield inc
  }

  implicit class EnhancedJson(val json: Json) extends AnyVal {
    def workflowId: ErrorOr[WorkflowId] = {
      val workflowIdOpt = for {
        o <- json.asObject
        id <- o.kleisli("id")
        s <- id.asString
      } yield WorkflowId.fromString(s)

      workflowIdOpt match {
        case Some(workflowId) => workflowId.validNel
        case None => s"did not find workflow id in ${json.printWith(Printer.spaces2).elided(100)}".invalidNel
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
  def updateLabels(json: Json, databaseLabels: Map[WorkflowId, Map[String, String]]): ErrorOr[Json] = {

    def updateLabelsInCalls(callObject: JsonObject, subworkflowJson: Json): Option[ErrorOr[Json]] = {
      // If the call contains a subWorkflowMetadata key, return a copy of the call with
      // its subworkflowMetadata updated via a recursive call to `doUpdateWorkflow`.
      Option(
        doUpdateWorkflow(subworkflowJson) map { updatedSubworkflow =>
          Json.fromJsonObject(callObject.add(subWorkflowMetadataKey, updatedSubworkflow))
        }
      )
    }

    def doUpdateWorkflow(workflowJson: Json): ErrorOr[Json] = {
      for {
        id <- workflowJson.workflowId
        workflowWithUpdatedCalls <- updateWorkflowCallsJson(workflowJson, updateLabelsInCalls)
        json = databaseLabels.get(id) match {
          case None => workflowWithUpdatedCalls
          case Some(labels) =>
            val labelsJson: Json = Json.fromFields(labels safeMapValues Json.fromString)
            workflowWithUpdatedCalls deepMerge Json.fromFields(List(("labels", labelsJson)))
        }
      } yield json
    }

    doUpdateWorkflow(workflowJson = json)
  }

  private def callsObject(workflowObject: JsonObject): ErrorOr[Option[JsonObject]] = {
    workflowObject("calls") match {
      case None => None.validNel
      case Some(callsJson) =>
        callsJson.asObject match {
          case None => s"'calls' member unexpectedly not a JSON object: $callsJson".invalidNel
          case Some(calls) => Option(calls).validNel
        }
    }
  }

  private def callsObject(workflowJson: Json): ErrorOr[Option[JsonObject]] = {
    workflowJson.asObject match {
      case None => s"Workflow JSON unexpectedly not an object: $workflowJson".invalidNel
      case Some(obj) => callsObject(obj)
    }
  }

  /**
    * Update workflow json, replacing "subWorkflowMetadata" elements of root workflow's "calls" object by "subWorkflowId"
    *
    * @param workflowJson json blob
    * @return updated json
    */
  def unexpandSubworkflows(workflowJson: Json): ErrorOr[Json] = {
    def workflowObject(workflowJson: Json): ErrorOr[JsonObject] =
      workflowJson.asObject.toErrorOr(s"Workflow JSON unexpectedly not a JSON object: $workflowJson")

    def workflowObjectWithUpdatedCalls(workflowObject: JsonObject, callsObject: JsonObject): ErrorOr[Json] = {
      def unexpandSubworkflow(callEntry: JsonObject, sub: Json): ErrorOr[Json] = {
        for {
          subObj <- sub.asObject.toErrorOr(s"subWorkflowMetadata unexpectedly not an object: $sub")
          id <- subObj("id").toErrorOr(s"subWorkflowMetadata unexpectedly missing 'id' field: $subObj")
          _ <- id.asString.toErrorOr(s"subworkflow 'id' field unexpectedly not a string: $id")
          updated = callEntry.remove(subWorkflowMetadataKey).add(subWorkflowIdKey, id)
        } yield Json.fromJsonObject(updated)
      }

      def possiblyUnexpandCallEntry(entry: Json): ErrorOr[Json] = {
        entry.asObject match {
          case None => s"call entry unexpectedly not an object: $entry".invalidNel
          case Some(entryObject) =>
            entryObject(subWorkflowMetadataKey) match {
              case None => entry.validNel // no subworkflow metadata is fine, return the entry unmodified.
              case Some(sub) => unexpandSubworkflow(entryObject, sub)
            }
        }
      }

      val updatedCallsObject = callsObject.traverse { json =>
        val updatedCallArray: ErrorOr[Vector[Json]] = json.asArray match {
          case Some(array) => array traverse possiblyUnexpandCallEntry
          case None => s"value unexpectedly not an array: $json".invalidNel
        }
        updatedCallArray map Json.fromValues
      }

      for {
        updatedCallsJson <- updatedCallsObject map Json.fromJsonObject
      } yield Json.fromJsonObject(workflowObject.add("calls", updatedCallsJson))
    }

    for {
      workflowObject <- workflowObject(workflowJson)
      callsObject <- callsObject(workflowObject)
      // `callsObject` is of type `Option[JsonObject]`, so the `map` has a return type of `Option[ErrorOr[Json]]`.
      // The `getOrElse` will return an expression of `ErrorOr[Json]`.
      updated <- callsObject map { workflowObjectWithUpdatedCalls(workflowObject, _) } getOrElse workflowJson.validNel
    } yield updated
  }

  private def updateWorkflowCallsJson(workflowJson: Json, updateCallsFunc: (JsonObject, Json) => Option[ErrorOr[Json]]): ErrorOr[Json] = {
    val workflowWithUpdatedCalls: ErrorOr[Json] = extractJsonObjectByKey(workflowJson, "calls") match {
      // If there were no calls just return the workflow JSON unmodified.
      case None => workflowJson.validNel
      case Some((_, calls)) =>
        val updatedCalls: ErrorOr[JsonObject] = calls.traverse {
          // The Json (a JSON array, really) corresponding to the array of call objects for a call name.
          callValue: Json =>
            // The object above converted to a Vector[Json].
            val callArray: Vector[Json] = callValue.asArray.toVector.flatten

            val updatedCallArray: Vector[Option[ErrorOr[Json]]] = callArray map { callJson =>
              // If there is no subworkflow object this will be None.
              val callAndSubworkflowObjects: Option[(JsonObject, Json)] = extractJsonByKey(callJson, subWorkflowMetadataKey)

              callAndSubworkflowObjects match {
                case None => Option(callJson.validNel)
                case Some((callObject, subworkflowJson)) => updateCallsFunc(callObject, subworkflowJson)
              }
            }
            val includedCallsOnly = updatedCallArray.flatten
            (includedCallsOnly.sequence[ErrorOr, Json]: ErrorOr[Vector[Json]]) map Json.fromValues
        }

        for {
          calls <- updatedCalls
          obj <- workflowJson.asObject.toErrorOr(s"unexpectedly not a JSON object: $workflowJson")
          updatedWorkflow = Json.fromJsonObject(obj.add("calls", Json.fromJsonObject(calls)))
        } yield updatedWorkflow
    }
    workflowWithUpdatedCalls
  }

  def extractSubWorkflowMetadata(subworkflowId: String, workflowJson: Json): ErrorOr[Option[Json]] = {
    extractJsonObjectByKey(workflowJson, "calls") match {
      case None => None.validNel
      case Some((_, calls)) =>
        calls.toMap.map {
          // The Json (a JSON array, really) corresponding to the array of call objects for a call name.
          case (_, callValue: Json) =>
            // The object above converted to a Vector[Json].
            val callArray: Vector[Json] = callValue.asArray.toVector.flatten

            // actually while this vector will contain one element per each object from `calls` array, only at most one
            // of those elements will be `Valid[Some[Json]]' (if subworkflow was found), while others should be `Valid[None]`
            val subworkflowsArrayFromCalls: Vector[ErrorOr[Option[Json]]] = callArray map { callJson =>
              // If there is no subworkflow object this will be None.
              val callAndSubworkflowObjects: Option[(JsonObject, Json)] = extractJsonByKey(callJson, subWorkflowMetadataKey)

              callAndSubworkflowObjects match {
                case None => None.validNel
                case Some((_, subworkflowJson)) =>
                  subworkflowJson.workflowId match {
                    case Valid(currentSubworkflowId) if currentSubworkflowId.toString == subworkflowId => Option(subworkflowJson).validNel
                    case Valid(_) => extractSubWorkflowMetadata(subworkflowId, subworkflowJson)
                    case err@Invalid(_) => err
                  }
              }
            }
            (subworkflowsArrayFromCalls.sequence[ErrorOr, Option[Json]]: ErrorOr[Vector[Option[Json]]]) map (_.flatten.headOption)
        }
          .toVector
          .sequence[ErrorOr, Option[Json]]
          .map(_.flatten.headOption)
    }
  }

}
