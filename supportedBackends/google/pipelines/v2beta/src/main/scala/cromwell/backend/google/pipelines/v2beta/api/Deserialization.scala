package cromwell.backend.google.pipelines.v2beta.api

import java.io.ByteArrayInputStream
import java.nio.charset.StandardCharsets
import java.util.{ArrayList => JArrayList, Map => JMap}
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.client.json.{GenericJson, JsonObjectParser}
import com.google.api.services.lifesciences.v2beta.model._
import common.validation.ErrorOr._
import common.validation.Validation._
import mouse.all._

import scala.collection.JavaConverters._
import scala.reflect.ClassTag
import scala.util.Try

/**
  * This bundles up some code to work around the fact that Operation is not deserialized
  * completely from the JSON HTTP response.
  * For instance, the metadata field is a Map[String, Object] even though it represents "Metadata"
  * for which there's an existing class.
  * This class provides implicit functions to deserialize those map to their proper type.
  */
private [api] object Deserialization {

  private val jsonFactory = new JacksonFactory
  private val jsonParser = new JsonObjectParser.Builder(jsonFactory).build

  implicit class OperationDeserialization(val operation: Operation) extends AnyVal {
    /**
      * Deserializes the events to com.google.api.services.genomics.v2beta.model.Event
      */
    def events: ErrorOr[List[Event]] = {
      val eventsErrorOrOption = for {
        eventsMap <- metadata.get("events")
        eventsErrorOr <- Option(eventsMap
          .asInstanceOf[JArrayList[JMap[String, Object]]]
          .asScala
          .toList
          .traverse[ErrorOr, Event](deserializeTo[Event](_).toErrorOr)
        )
      } yield eventsErrorOr
      eventsErrorOrOption.getOrElse(Nil.validNel)
    }

    /**
      * Deserializes the pipeline to com.google.api.services.genomics.v2beta.model.Pipeline
      */
    def pipeline: Option[Try[Pipeline]] = {
      metadata
        .get("pipeline")
        .map(_.asInstanceOf[JMap[String, Object]] |> deserializeTo[Pipeline])
    }

    // If there's a WorkerAssignedEvent it means a VM was created - which we consider as the job started
    // Note that the VM might still be booting
    def hasStarted = events.toOption.exists(_.exists(_.getWorkerAssigned != null))

    def metadata: Map[String, AnyRef] = {
      val metadataOption = for {
        operationValue <- Option(operation)
        metadataJava <- Option(operationValue.getMetadata)
      } yield metadataJava.asScala.toMap
      metadataOption.getOrElse(Map.empty)
    }
  }

  /**
    * Deserializes a java.util.Map[String, Object] to an instance of T
    */
  private [api] def deserializeTo[T <: GenericJson](attributes: JMap[String, Object])(implicit tag: ClassTag[T]): Try[T] = Try {
    val jsonString = jsonFactory.toString(attributes)
    jsonParser.parseAndClose[T](new ByteArrayInputStream(jsonString.getBytes), StandardCharsets.UTF_8, tag.runtimeClass.asInstanceOf[Class[T]])
  }
}
