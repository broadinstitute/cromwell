package cromwell.backend.impl.htcondor.caching.provider.mongodb.model

import spray.json

/**
  * Wrapper over the byte array that is stored in db
  *
  * @param byteArray Serialized data
  */
case class KryoSerializedObject(byteArray: Array[Byte])

/**
  * SucceededResponse to be stored in MongoDB.
  *
  * @param hash Calculated hash for the Job.
  * @param succeededResponse Serialized succeeded response.
  */
case class MongoCachedExecutionResult(hash: String, succeededResponse: KryoSerializedObject)

object MongoCachedExecutionResultProtocol extends json.DefaultJsonProtocol {
  implicit val kryoSerializedObject = jsonFormat1(KryoSerializedObject)
  implicit val cachedExecutionResultProtocol = jsonFormat2(MongoCachedExecutionResult)
}

