package cromwell.backend.google.pipelines.v2alpha1.api

import java.util

import com.google.api.services.genomics.model.WorkerAssignedEvent
import com.google.api.services.genomics.v2alpha1.model.{ContainerStartedEvent, Operation}
import cromwell.backend.google.pipelines.v2alpha1.api.Deserialization._
import org.scalatest.{FlatSpec, Matchers}

import scala.collection.JavaConverters._

class DeserializationSpec extends FlatSpec with Matchers {
  behavior of "Deserialization"

  it should "deserialize events from operation metadata" in {
    val operation = new Operation()
    val metadataMap = Map[String, AnyRef](
      "events" -> new util.ArrayList(
        List[java.util.Map[String, Object]](
          Map[String, AnyRef](
            "description" -> "event 1 description",
            "timestamp" -> "2018-04-20T14:38:25+00:00",
            "details" -> Map[String, Object](
              "@type" -> "WorkerAssignedEvent",
              "zone" -> "event 1 Zone",
              "instance" -> "event 1 Instance"
            ).asJava
          ).asJava,
          Map[String, AnyRef](
            "description" -> "event 2 description",
            "timestamp" -> "2018-04-20T14:39:25+00:00",
            "details" -> Map[String, Object](
              "@type" -> "ContainerStartedEvent",
              "actionId" -> Integer.valueOf(18),
              "ipAddress" -> "86.127.54.8",
              "portMappings" -> Map(
                8000 -> 8008
              ).asJava
            ).asJava
          ).asJava
        ).asJava
      )
    ).asJava

    operation.setMetadata(metadataMap)
    val deserializedEvents = operation.events.valueOr(errors => fail(errors.toList.mkString(", ")))

    val event1 = deserializedEvents.head
    event1.getDescription shouldBe "event 1 description"
    event1.getTimestamp shouldBe "2018-04-20T14:38:25+00:00"
    // Event1 details are of type WorkerAssignedEvent, so it should not be defined for something else 
    event1.details[ContainerStartedEvent] should not be defined

    val event1Details = event1.details[WorkerAssignedEvent].map(_.get)
    event1Details shouldBe defined
    event1Details.get.getInstance shouldBe "event 1 Instance"
    event1Details.get.getZone shouldBe "event 1 Zone"

    val event2 = deserializedEvents(1)
    event2.getDescription shouldBe "event 2 description"
    event2.getTimestamp shouldBe "2018-04-20T14:39:25+00:00"

    val event2Details = event2.details[ContainerStartedEvent].map(_.get)
    event2Details shouldBe defined
    event2Details.get.getActionId shouldBe 18
    event2Details.get.getIpAddress shouldBe "86.127.54.8"
    event2Details.get.getPortMappings.size() shouldBe 1
    event2Details.get.getPortMappings.get(8000) shouldBe 8008
  }

  it should "deserialize pipeline from operation metadata" in {
    val operation = new Operation()

    val metadataMap = Map[String, AnyRef](
      "pipeline" -> Map[String, AnyRef](
        "actions" -> List[java.util.Map[String, Object]](
          Map[String, Object](
            "name" -> "actionName",
            "imageUri" -> "ubuntu:latest",
            "commands" -> List[String]("echo", "hello").asJava
          ).asJava
        ).asJava,
        "resources" -> Map(
          "projectId" -> "project",
          "virtualMachine" -> Map(
            "machineType" -> "custom-1-1024",
            "preemptible" -> false
          ).asJava
        ).asJava
      ).asJava
    ).asJava

    operation.setMetadata(metadataMap)
    val deserializedPipeline = operation.pipeline.get.get
    val action = deserializedPipeline.getActions.get(0)
    action.getCommands.asScala shouldBe List("echo", "hello")
    action.getImageUri shouldBe "ubuntu:latest"
    action.getName shouldBe "actionName"
    deserializedPipeline.getResources.getProjectId shouldBe "project"
    val virtualMachine = deserializedPipeline.getResources.getVirtualMachine
    virtualMachine.getMachineType shouldBe "custom-1-1024"
    virtualMachine.getPreemptible shouldBe false
  }

  // https://github.com/broadinstitute/cromwell/issues/4772
  it should "deserialize pipeline from operation metadata without preemptible" in {
    val operation = new Operation()

    val metadataMap = Map[String, AnyRef](
      "pipeline" -> Map[String, AnyRef](
        "actions" -> List[java.util.Map[String, Object]](
          Map[String, Object](
            "name" -> "actionName",
            "imageUri" -> "ubuntu:latest",
            "commands" -> List[String]("echo", "hello").asJava
          ).asJava
        ).asJava,
        "resources" -> Map(
          "projectId" -> "project",
          "virtualMachine" -> Map(
            "machineType" -> "custom-1-1024",
          ).asJava
        ).asJava
      ).asJava
    ).asJava

    operation.setMetadata(metadataMap)
    val deserializedPipeline = operation.pipeline.get.get
    val action = deserializedPipeline.getActions.get(0)
    action.getCommands.asScala shouldBe List("echo", "hello")
    action.getImageUri shouldBe "ubuntu:latest"
    action.getName shouldBe "actionName"
    deserializedPipeline.getResources.getProjectId shouldBe "project"
    val virtualMachine = deserializedPipeline.getResources.getVirtualMachine
    virtualMachine.getMachineType shouldBe "custom-1-1024"
    virtualMachine.getPreemptible shouldBe null
  }

  it should "be able to say if the operation has started" in {
    val operation = new Operation()
    
    def makeMetadata(details: Map[String, Object]) = Map[String, AnyRef](
      "events" -> new util.ArrayList(
        List[java.util.Map[String, Object]](
          Map[String, AnyRef](
            "description" -> "event 1 description",
            "timestamp" -> "2018-04-20T14:38:25+00:00",
            "details" -> details.asJava
          ).asJava
        ).asJava
      )
    ).asJava
    
    val metadataMapStarted = makeMetadata(Map[String, Object](
      "@type" -> "WorkerAssignedEvent",
      "zone" -> "event 1 Zone",
      "instance" -> "event 1 Instance"
    ))
    val metadataMapNotStarted = makeMetadata(Map.empty)
    val metadataMapNotStarted2 = makeMetadata(Map[String, Object](
      "@type" -> "ContainerStartedEvent"
    ))
    
    operation.setMetadata(metadataMapStarted)
    operation.hasStarted shouldBe true
    operation.setMetadata(metadataMapNotStarted)
    operation.hasStarted shouldBe false
    operation.setMetadata(metadataMapNotStarted2)
    operation.hasStarted shouldBe false
  }
  
  it should "deserialize big decimals correctly" in {
    val valueMap = Map[String, Object](
      "integerValue" -> BigDecimal(5),
      "doubleValue" -> BigDecimal.decimal(6D),
      "floatValue" -> BigDecimal.decimal(7F),
      "longValue" -> BigDecimal.decimal(8L)
    ).asJava
    
    val deserialized = Deserialization.deserializeTo[DeserializationTestClass](valueMap)
    deserialized.isSuccess shouldBe true
    val deserializedSuccess = deserialized.get
    deserializedSuccess.integerValue shouldBe 5
    deserializedSuccess.doubleValue shouldBe 6D
    deserializedSuccess.floatValue shouldBe 7F
    deserializedSuccess.longValue shouldBe 8L
  }
  
}
