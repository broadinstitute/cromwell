package cromwell.backend.impl.htcondor.caching.provider.mongodb.serialization

import java.io.{ByteArrayInputStream, ByteArrayOutputStream}

import com.esotericsoftware.kryo.io.{Input, Output}
import com.twitter.chill.ScalaKryoInstantiator

/**
  * This mixin provides access to `serialize` and `deserialize` methods that use Kryo underneath to
  * perform the actual conversion to / from bytes
  */
trait KryoSerDe extends SerDe {

  /**
    *
    * @param data Any Scala / Java object that needs to be serialized
    * @return A serialized byte array that can be transported across Network moved around with
    */
  override def serialize[A <: AnyRef](data: A): Array[Byte] = {
    try {
      val instantiator = new ScalaKryoInstantiator
      instantiator.setRegistrationRequired(false) // This makes it unnecessary to register all classes
      val kryo = instantiator.newKryo()
      val buffer = new ByteArrayOutputStream()
      val output = new Output(buffer)
      kryo.writeObject(output, data)
      output.close()
      buffer.toByteArray
    } catch {
      case exception: Exception => throw SerDeException("Failed to serialize data.", exception)
    }
  }

  /**
    *
    * @param byteArray Kryo serialized data. Expects only results of `writeObject`
    * @param toClass Class to which the `byteArray` will be deserialized
    * @return The deserialized object as an instance of A
    */
  override def deserialize[A <: AnyRef](byteArray: Array[Byte], toClass: Class[A]): A = {
    try {
      val instantiator = new ScalaKryoInstantiator
      instantiator.setRegistrationRequired(false)
      val kryo = instantiator.newKryo()
      val buffer = new ByteArrayInputStream(byteArray)
      val input = new Input(buffer)
      val result = kryo.readObject(input, toClass)
      input.close()
      result
    } catch {
      case exception: Exception => throw SerDeException("Failed to deserialize data.", exception)
    }
  }

}
