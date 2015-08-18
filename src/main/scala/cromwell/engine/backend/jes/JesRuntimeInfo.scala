package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.{DockerExecutor, Resources}
import cromwell.binding.{Call, RuntimeAttributes}

import scala.collection.JavaConverters._

object JesRuntimeInfo {
  def apply(commandLine: String, call: Call): JesRuntimeInfo = {
    val dockerImage = call.docker.getOrElse(throw new IllegalArgumentException("FIXME: WHAT TO DO?"))

    new JesRuntimeInfo(buildResources(call), buildDockerExecutor(commandLine, dockerImage))
  }

  def buildDockerExecutor(commandLine: String, dockerImage: String): DockerExecutor = {
    val docker = new DockerExecutor()
    docker.setImage(dockerImage).setCmd(commandLine)
  }

  def buildResources(call: Call): Resources = {
    val runtimeAttributes: RuntimeAttributes = call.task.runtimeAttributes

    new Resources()
      .setRamGb(runtimeAttributes.memoryGB)
      .setCpu(runtimeAttributes.cpu)
      .setPreemptible(runtimeAttributes.preemptible)
      .setZones(runtimeAttributes.defaultZones.asJava)
      .setDisks(runtimeAttributes.defaultDisks.asJava)
  }
}

case class JesRuntimeInfo(resources: Resources, docker: DockerExecutor)
