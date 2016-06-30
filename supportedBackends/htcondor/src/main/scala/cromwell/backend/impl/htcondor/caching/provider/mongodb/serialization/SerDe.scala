package cromwell.backend.impl.htcondor.caching.provider.mongodb.serialization

/**
  * A vanilla Serialization / Deserialization interface
  */
trait SerDe {

  case class SerDeException(message: String, error: Throwable) extends IllegalStateException(message, error)

  /**
    * Serialize the given data
    * @param data Any Scala / Java object that needs to be serialized
    * @return A serialized byte array that can be transported across Network moved around with
    */
  def serialize[A <: AnyRef](data: A): Array[Byte]

  /**
    * Deserializes an array of Bytes to the given class
    * @param byteArray Kryo serialized data. Expects only results of `writeObject`
    * @param toClass Class to which the `byteArray` will be deserialized
    * @return The deserialized object as an instance of A
    */
  def deserialize[A <: AnyRef](byteArray: Array[Byte], toClass: Class[A]): A

}
