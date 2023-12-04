package cromwell.engine.workflow.lifecycle.execution.job.preparation

import cromwell.services.keyvalue.KeyValueServiceActor.{KvResponse, ScopedKey}

/**
  * Handles the determination of when we know key lookups are successful.
  */
sealed private trait KeyValueLookups

final private[preparation] case class PartialKeyValueLookups(responses: Map[ScopedKey, KvResponse],
                                                             awaiting: Seq[ScopedKey]
) {
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

final private case class KeyValueLookupResults(values: Map[ScopedKey, KvResponse]) {
  def unscoped: Map[String, KvResponse] = values map { case (k, v) => k.key -> v }
}
