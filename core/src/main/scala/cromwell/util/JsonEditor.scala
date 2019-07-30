package cromwell.util

import cats.data.NonEmptyList
import common.collections.EnhancedCollections._
import io.circe.{ACursor, FailedCursor}
import io.circe.Json
import mouse.all._
import cats.syntax.traverse._
import cats.syntax.functor._
import cats.instances.list._
import cats.data.State

object JsonEditor {

  def includeExcludeJson(json: Json, includeKeys: Option[NonEmptyList[String]], excludeKeys: Option[NonEmptyList[String]]): Json = {
    (includeKeys, excludeKeys) match {
      // Take includes, and then remove excludes
      case (Some(includeKeys), Some(excludeKeys)) => includeJson(json, includeKeys) |> (excludeJson(_, excludeKeys))
      case (None, Some(excludeKeys)) => excludeJson(json, excludeKeys)
      case (Some(includeKeys), None) => includeJson(json, includeKeys)
      case _ => json
    }
  }

  private def excludeKeys(keys: NonEmptyList[String]) : State[ACursor, Unit] = State.modify{cursor =>
    cursor.keys match {
      case Some(levelKeys) =>
        val modifications = levelKeys.toList.traverse{key =>
          val delete = keys.foldLeft(false)(_ || key.startsWith(_))
          State.modify[ACursor] { cursor =>
            if (delete)
              //moves cursor back to parent
              cursor.downField(key).delete
            else {
              val newCursor = cursor.downField(key)
              val eval = excludeKeys(keys).run(newCursor)
              //ignoring unit output
              val (output,_) = eval.value
              //we're in the field, have to go back to the parent for the next field to evaluate
              output.up
            }
          }
        }.void
        val (newCursor, _) = modifications.run(cursor).value
        newCursor
      case _ =>
        val arrayFirstElement = cursor.downArray
        arrayFirstElement match {
          case _: FailedCursor => cursor
          case _ =>
            val (newCursor,_) = excludeKeysArray(keys).run(arrayFirstElement).value
            //go back to our parent as above
            newCursor.up
        }

    }
  }
  private def excludeKeysArray(keys: NonEmptyList[String]) : State[ACursor, Unit] = State.modify {
    cursor =>
      val (newCursor,_) = excludeKeys(keys).run(cursor).value
      newCursor.right match {
        case _:FailedCursor => newCursor
        case rightCursor =>
          val (nextState,_) = excludeKeysArray(keys).run(rightCursor).value
          nextState
      }
  }

  /**
    * @param keys list of keys to match against field names using "startsWith."
    * @return A state transition that will tell whether or not to keep this path in order to keep a nested value.
    */
  private def includeKeys(keys: NonEmptyList[String]) : State[ACursor, Boolean] = State.apply[ACursor, Boolean]{cursor =>
    cursor.keys match {
      case Some(levelKeys) =>
        val modifications: State[ACursor, List[Boolean]] = levelKeys.toList.traverse { key =>
          //detect whether any of the json keys match the argument keys
          val keep = keys.foldLeft(false)((acc, i) => acc || key.startsWith(i))
          State.apply[ACursor, Boolean] { cursor =>
            if (keep) {
              (cursor, true)
            } else {
              val newCursor = cursor.downField(key)

              //before deleting, we want to see if any children should be kept.  The boolean output will tell us that
              val (output, shouldKeep) = includeKeys(keys).run(newCursor).value

              if (shouldKeep) {
                (output.up, true)
              } else {
                //no need to keep this node on account of a child needing it.
                //delete this node and indicate to our own parent that we don't need it.
                (newCursor.delete, false)
              }
            }
          }
        }
        val (newCursor, keeps) = modifications.run(cursor).value
        (newCursor, keeps.foldLeft(false)(_ || _))
      case None =>
        val arrayFirstElement = cursor.downArray
        arrayFirstElement match {
          case _:FailedCursor => (cursor, false)
          case _ =>
            val (nextState,keep) = includeKeysArray(keys).run(arrayFirstElement).value
            (nextState.up, keep)
        }
    }
  }

  private def includeKeysArray(keys: NonEmptyList[String]) : State[ACursor, Boolean] = State.apply[ACursor, Boolean] {
    cursor =>
      val (newCursor,keep) = includeKeys(keys).run(cursor).value
      newCursor.right match {
        case _:FailedCursor => (newCursor, keep)
        case rightCursor =>
          val (nextState,keep2) = includeKeysArray(keys).run(rightCursor).value
          (nextState, keep || keep2)
      }
  }

  private def modifyJson[A](json: Json, keys: NonEmptyList[String],  function: NonEmptyList[String] => State[ACursor, A]) = {
    val cursor = json.hcursor
    val mod: State[ACursor, A] = function(keys)
    val (newCursor, _) = mod.run(cursor).value
    //taking the liberty of assuming the document still has something in it.  Arguable, might warrant further consideration.
    newCursor.top.get
  }

  def includeJson(json: Json, keys: NonEmptyList[String]) = modifyJson(json, keys, includeKeys)

  def excludeJson(json: Json, keys: NonEmptyList[String]) = modifyJson(json, keys, excludeKeys)

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
