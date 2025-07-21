package cromwell.services.keyvalue.impl

import akka.actor.{ActorRef, Props}
import cromwell.services.keyvalue.{KeyValueServiceActor, KeyValueWriteActor}

import scala.concurrent.duration.FiniteDuration

object SqlKeyValueWriteActor {
  def props(threshold: Int, serviceRegistryActor: ActorRef, flushRate: FiniteDuration, batchSize: Int) =
    Props(new SqlKeyValueWriteActor(threshold, serviceRegistryActor, flushRate, batchSize))
}

class SqlKeyValueWriteActor(override val threshold: Int,
                            serviceRegistryActor: ActorRef,
                            flushRate: FiniteDuration,
                            batchSize: Int
) extends KeyValueWriteActor(serviceRegistryActor, flushRate, batchSize)
    with BackendKeyValueDatabaseAccess {
  implicit private val system = context.system

  override def processPut(puts: Vector[KeyValueServiceActor.KvPut]) = {
    val pairs = puts.map { put =>
      (put.pair.key.workflowId, put.pair.key.jobKey, put.pair.key.key, put.pair.value)
    }
    updateBackendKeyValuePairs(pairs)
  }
}
