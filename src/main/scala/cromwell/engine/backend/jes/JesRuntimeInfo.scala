package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.{DockerExecutor, Resources}
import cromwell.binding.{Call, RuntimeAttributes}

import scala.collection.JavaConverters._

object JesRuntimeInfo {
  def apply(commandLine: String, call: Call): JesRuntimeInfo = {
    /*
     It should be impossible for docker to be None here. Enforcing that w/ ADTs seemed more trouble than
      it was worth. If you're ever debugging a NoSuchElementException which leads you here, that means
      the more trouble than worth calculation was incorrect and we should have separate RuntimeAttributes for
      docker and no docker cases
    */
    val dockerImage = call.docker.get

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
