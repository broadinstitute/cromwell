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
import io.circe.{Json, JsonObject, Printer}

object JsonEditor {

  private val subWorkflowMetadataKey = "subWorkflowMetadata"
  private val subWorkflowIdKey = "subWorkflowId"

  def includeExcludeJson(json: Json, includeKeys: Option[NonEmptyList[String]], excludeKeys: Option[NonEmptyList[String]]): ErrorOr[Json] =
    (includeKeys, excludeKeys) match {
      // Take includes, and then remove excludes
      case (Some(includeKeys), Some(excludeKeys)) => includeJson(json, includeKeys) flatMap (excludeJson(_, excludeKeys))
      case (None, Some(excludeKeys)) => excludeJson(json, excludeKeys)
      case (Some(includeKeys), None) => includeJson(json, includeKeys)
      case _ => json.validNel
    }

  final case class Filter(components: NonEmptyList[String])

  /** A `FilterGroup` represents all the `Filter`s for include xor exclude. The argument is intentionally not a NEL
    * since `FilterGroup`s should start out non-empty but may become empty after `remove`s. e.g. excludeKey 'id' in
    * a workflow `FilterGroup`. */
  final case class FilterGroup(filters: List[Filter]) {
    def withoutFilters(h: String, t: String*): FilterGroup = {
      val componentSet = (h :: t.toList).toSet
      this.copy(filters filterNot { f => componentSet.contains(f.components.head) })
    }

    private def filtersMatchingKey(objectKey: String): List[Filter] = filters.filter(_.components.head == objectKey)

    def applyExcludes(jsonObject: JsonObject): JsonObject = {
      // If at the top level there are no matches at all, either because there is no object to be drilled into or
      // there is no object key matching the current filter component, keep the original JSON.
      // If the key has matched completely, remove the key/value from the object and return the modified object.
      // Otherwise recursively descend into the object attempting to continue to match.
      def descend(json: Json, remainingComponents: NonEmptyList[String]): Json = {
        json.asObject match {
          case None => json
          case Some(jsonObject) =>
            val component = remainingComponents.head
            jsonObject(component) match {
              case None => json
              case Some(inner) =>
                val updatedObject = remainingComponents.tail match {
                  case Nil =>
                    jsonObject.remove(component)
                  case h :: t =>
                    val descentResult = descend(inner, NonEmptyList.of(h, t: _*))
                    val innerObjectBecameEmpty = for {
                      innerObject <- inner.asObject
                      descentResultObject <- descentResult.asObject
                    } yield innerObject.nonEmpty && descentResultObject.isEmpty

                    // Special case: if the JSON object has become empty as a result of editing, delete it.
                    if (innerObjectBecameEmpty.contains(true))
                      jsonObject.remove(component)
                    else
                      jsonObject.add(component, descentResult)
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
            val filteredJson = fs.foldLeft(json) { case (j, f) =>
              f.components.tail match {
                // We know all of these filters had "tails". i.e. they did not match completely.
                case Nil => throw new RuntimeException("Programmer error: filter components tail should not be empty")
                case c :: cs => descend(j, NonEmptyList.of(c, cs: _*))
              }
            }
            List((key, filteredJson))
        }
      }
      JsonObject.fromIterable(updatedKeyValues)
    }

    def applyIncludes(jsonObject: JsonObject): JsonObject = {
      // If the JSON is not an object or there is no key in the object matching the head of the remaining components,
      // the search has failed so return `None`. If there is a matching key and there are no more components,
      // the search has succeeded so return a result. Otherwise continue the search deeper into the inner JSON value.
      def descend(json: Json, remainingComponents: NonEmptyList[String]): Option[Json] = {
        json.asObject match {
          case None => None
          case Some(jsonObject) =>
            val component = remainingComponents.head
            jsonObject(component) match {
              case None => None
              case Some(inner) =>
                remainingComponents.tail match {
                  case Nil =>
                    // Pick out the matching key / value pair from the query object and return that as a
                    // JSON object with a single key / value.
                    Option(Json.fromFields(List((component, inner))))
                  case h :: t =>
                    val descentResult = descend(inner, NonEmptyList.of(h, t: _*))
                    // If there is a result, wrap it in a single field JSON object using the head of the filter
                    // components as the key. This rebuilds the object structure from the inside out as the
                    // recursion backs out.
                    descentResult map { j => Json.fromFields(List((h, j)))}
                }
            }
        }
      }

      val updatedKeyValues: List[(String, Json)] = jsonObject.toList flatMap { case (key, json) =>
        filtersMatchingKey(key) match {
          case Nil =>
            // No filters matched this object key, return an empty List.
            Nil
          case fs if fs.exists(_.components.tail == Nil) =>
            // If there is a filter that has no other components and thus matches the object key completely,
            // return the (key, value).
            List((key, json))
          case fs =>
            // The search continues. Fold an empty JsonObject through all the `fs` filters to include everything that matches.
            val emptyJson = Json.fromJsonObject(JsonObject.empty)
            val filteredJson = fs.foldLeft(emptyJson) { case (j, f) =>
              f.components.tail match {
                // We know all of these filters had "tails". i.e. they did not match completely.
                case Nil => throw new RuntimeException("Programmer error: filter components tail should not be empty")
                case c :: cs =>
                  descend(json, NonEmptyList.of(c, cs: _*)) match {
                    case None =>
                      // If the descent doesn't find anything return the object under construction unmodified.
                      j
                    case Some(descentResult) =>
                      // If the descent does find a value it should be merged to the object under construction.
                      // Merge is required to properly handle nested structure if there are multiple semi-overlapping
                      // include filter keys.
                      j deepMerge descentResult
                  }
              }
            }
            List((key, filteredJson))
        }
      }
      JsonObject.fromIterable(updatedKeyValues)
    }
  }

  object FilterGroup {
    def build(keys: NonEmptyList[String]): ErrorOr[FilterGroup] = {
      val filters: ErrorOr[List[Filter]] = keys.toList traverse { key => key.split(':').toList match {
        case Nil => s"Programmer error: string split resulting in empty array: $key".invalidNel
        case h :: t => Filter(NonEmptyList.of(h, t: _*)).validNel
      }}
      filters map FilterGroup.apply
    }
  }

  def excludeJson(json: Json, keys: NonEmptyList[String]): ErrorOr[Json] = FilterGroup.build(keys) flatMap applyExcludes(json)

  def includeJson(json: Json, keys: NonEmptyList[String]): ErrorOr[Json] = FilterGroup.build(keys) flatMap applyIncludes(json, forceWorkflowId = true)

  private def applyIncludes(workflowJson: Json, forceWorkflowId: Boolean = false)(filterGroup: FilterGroup): ErrorOr[Json] = {
    // Returns an empty JsonObject if all calls have been filtered out or if there were never any calls to begin with.
    def filterCalls(workflowObject: JsonObject): ErrorOr[JsonObject] = {
      def buildCallObject(unfilteredCallObject: JsonObject, shallowFilteredCallObject: JsonObject, filteredSubworkflow: JsonObject): JsonObject = {
        // Add in the subworkflow if it is nonempty after filtering, add in `shardIndex` and `attempt` if anything is nonempty.
        val transforms: List[JsonObject => JsonObject] = List(
          obj =>
            if (filteredSubworkflow.nonEmpty)
              obj.add(subWorkflowMetadataKey, Json.fromJsonObject(filteredSubworkflow)) else obj,
          obj =>
            if (shallowFilteredCallObject.nonEmpty || filteredSubworkflow.nonEmpty)
              obj.add("shardIndex", unfilteredCallObject("shardIndex").get).add("attempt", unfilteredCallObject("attempt").get) else obj
        )
        Function.chain(transforms)(shallowFilteredCallObject)
      }

      def filterCallEntry(json: Json): ErrorOr[JsonObject] = {
        for {
          callObject <- json.asObject.map(_.validNel).getOrElse(s"Call entry unexpectedly not an object: $json".invalidNel)
          shallowFilteredCallObject = filterGroup.applyIncludes(callObject)
          subworkflow = callObject(subWorkflowMetadataKey).getOrElse(Json.fromJsonObject(JsonObject.empty))
          shallowFilteredSubworkflow <- applyIncludes(subworkflow)(filterGroup)
          shallowFilteredSubworkflowObject <- shallowFilteredSubworkflow.asObject.map(_.validNel).getOrElse(s"subworkflow unexpectedly not an object: $subworkflow ".invalidNel)
          builtCallObject = buildCallObject(callObject, shallowFilteredCallObject, shallowFilteredSubworkflowObject)
        } yield builtCallObject
      }

      def filterCallsArray(json: Json): ErrorOr[Json] = {
        json.asArray match {
          case None => s"Calls array unexpectedly not an array: $json".invalidNel
          case Some(array) =>
            val filteredCalls: ErrorOr[Vector[JsonObject]] = array traverse filterCallEntry
            // Remove empty call objects from the call array.
            val nonEmptyObjects: ErrorOr[Vector[JsonObject]] = filteredCalls map { _.filter(_.nonEmpty) }
            nonEmptyObjects map { objs => Json.fromValues(objs map Json.fromJsonObject) }
        }
      }

      val callsAsObject: ErrorOr[JsonObject] = workflowObject("calls") match {
        case None => JsonObject.empty.validNel
        case Some(cs) =>
         cs.asObject map { _.validNel } getOrElse s"calls JSON unexpectedly not an object: $cs".invalidNel
      }

      for {
        callsObject <- callsAsObject
        // Empty calls filtered out of call arrays.
        emptyCallsRemoved <- callsObject traverse filterCallsArray
        // Empty call arrays removed from the containing calls object.
        emptyCallArraysRemoved = emptyCallsRemoved filter { case (_, json) => json.asArray.exists(_.nonEmpty) }
      } yield emptyCallArraysRemoved
    }

    def buildWorkflowObject(workflowObject: JsonObject,
                            filteredCalls: JsonObject,
                            filteredWorkflow: JsonObject,
                            forceWorkflowId: Boolean): JsonObject = {

      // Add in `calls` if they are nonempty, add in `id` if anything is nonempty or the workflow ID is being forced
      // in as a result of being the query workflow ID.
      val transforms: List[JsonObject => JsonObject] = List(
        obj =>
          if (filteredCalls.nonEmpty)
            obj.add("calls", Json.fromJsonObject(filteredCalls)) else obj,
        obj =>
          if (filteredCalls.nonEmpty || filteredWorkflow.nonEmpty || forceWorkflowId)
            obj.add("id", workflowObject("id").get) else obj
      )
      Function.chain(transforms)(filteredWorkflow)
    }

    import common.validation.ErrorOr._
    for {
      workflowObject <- workflowJson.asObject.map(_.validNel).getOrElse(s"Workflow JSON unexpectedly not an object: $workflowJson".invalidNel)
      filteredCalls <- filterCalls(workflowObject)
      filteredWorkflow = filterGroup.applyIncludes(workflowObject)
      builtWorkflowObject = buildWorkflowObject(workflowObject, filteredCalls, filteredWorkflow, forceWorkflowId)
    } yield Json.fromJsonObject(builtWorkflowObject)
  }

  private def applyExcludes(workflowJson: Json)(filterGroup: FilterGroup): ErrorOr[Json] = {
    // Always include 'id' in workflows, 'shardIndex' and 'attempt' in calls. i.e. make sure these keys are excluded
    // from the excludes.
    val workflowExcludeFilters = filterGroup.withoutFilters("id")
    val callExcludeFilters = filterGroup.withoutFilters("shardIndex", "attempt")

    def filterCall(callJson: Json): ErrorOr[Json] = {
      def filterSubworkflowCall(parentWorkflowCall: JsonObject)(subworkflow: Json): ErrorOr[JsonObject] = {
        applyExcludes(subworkflow)(filterGroup) map { filteredSubworkflow => parentWorkflowCall.add(subWorkflowMetadataKey, filteredSubworkflow) }
      }

      for {
        callObject <- callJson.asObject.map(_.validNel).getOrElse(s"call JSON unexpectedly not an object: $callJson".invalidNel)
        // This applies filters based only on the call's keys, not looking for subworkflows.
        shallowFilteredCallObject = callExcludeFilters.applyExcludes(callObject)
        // If the call array element represents a subworkflow call this will recursively apply workflow and call filters as appropriate.
        deepFilteredCallObject <- shallowFilteredCallObject(subWorkflowMetadataKey) map filterSubworkflowCall(shallowFilteredCallObject) getOrElse shallowFilteredCallObject.validNel
      } yield Json.fromJsonObject(deepFilteredCallObject)
    }

    def filterCallsArray(callsArrayJson: Json): ErrorOr[Json] = {
      for {
        callsArray <- callsArrayJson.asArray.map(_.validNel).getOrElse(s"calls JSON unexpectedly not an array: $callsArrayJson".invalidNel)
        filteredCalls <- callsArray.traverse[ErrorOr, Json](filterCall)
      } yield Json.fromValues(filteredCalls)
    }

    def filterCallsObject(jsonObject: JsonObject): ErrorOr[Json] =
      (jsonObject traverse filterCallsArray) map Json.fromJsonObject

    for {
      workflowObject <- workflowJson.asObject.map(_.validNel).getOrElse(s"Workflow JSON unexpectedly not an object: $workflowJson".invalidNel)
      // This applies filters based only on the workflow's keys.
      shallowFilteredWorkflowObject = workflowExcludeFilters.applyExcludes(workflowObject)
      deepFilteredWorkflowObject <- shallowFilteredWorkflowObject("calls") match {
        case None => shallowFilteredWorkflowObject.validNel
        case Some(calls) => calls.asObject match {
          case None => s"calls JSON unexpectedly not an object: $calls".invalidNel
          case Some(callsObject) => filterCallsObject(callsObject) map {
            filteredCallsObject => shallowFilteredWorkflowObject.add("calls", filteredCallsObject)
          }
        }
      }
    } yield Json.fromJsonObject(deepFilteredWorkflowObject)
  }

  def outputs(json: Json): ErrorOr[Json] = {
    includeExcludeJson(json, Option(NonEmptyList.one("outputs")), Option(NonEmptyList.one("calls")))
  }

  def logs(workflowJson: Json): ErrorOr[Json] = {
    for {
      unexpanded <- unexpandSubworkflows(workflowJson)
      logs <- includeJson(unexpanded, NonEmptyList.of("stdout", "stderr", "backendLogs"))
    } yield logs
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
    * @return - optional tuple of workflow JsonObject and found element JsonObject
    */
  private def extractCallsKeyValue(workflowJson: Json): Option[(JsonObject, JsonObject)] =
    extractJsonByKey(workflowJson, "calls") flatMap {
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
    val workflowWithUpdatedCalls: ErrorOr[Json] = extractCallsKeyValue(workflowJson) match {
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
    extractCallsKeyValue(workflowJson) match {
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
