package cromwell.engine.workflow.lifecycle.execution.preparation

import cromwell.services.keyvalue.KeyValueServiceActor.{KvResponse, ScopedKey}

/**
  * Handles the determination of when we know key lookups are successful.
  */
private sealed trait KeyValueLookups

private[preparation] final case class PartialKeyValueLookups(responses: Map[ScopedKey, KvResponse], awaiting: Seq[ScopedKey]) {
  def withResponse(key: ScopedKey, response: KvResponse) = {
    val newResponses = responses + (key -> response)
    val newAwaiting = awaiting diff List(key)
    if (newAwaiting.isEmpty) {
      KeyValueLookupResults(newResponses)
    } else {
      PartialKeyValueLookups(newResponses, newAwaiting)
    }
  }
}

private final case class KeyValueLookupResults(values: Map[ScopedKey, KvResponse]) {
  def unscoped: Map[String, KvResponse] = values map { case (k, v) => k.key -> v }
}
