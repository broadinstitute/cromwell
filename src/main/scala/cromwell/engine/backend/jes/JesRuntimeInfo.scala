package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.{DockerExecutor, Resources}
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes

import scala.collection.JavaConverters._

sealed trait JesRuntimeInfo {
  def resources: Resources
  def docker: DockerExecutor
}

trait JesRuntimeInfoBuilder {
  def buildDockerExecutor(commandLine: String, dockerImage: String): DockerExecutor = {
    val docker = new DockerExecutor()
    docker.setImage(dockerImage).setCmd(commandLine)
  }

  def buildResources(runtimeAttributes: CromwellRuntimeAttributes): Resources = {
    new Resources()
      .setRamGb(runtimeAttributes.memoryGB)
      .setCpu(runtimeAttributes.cpu)
      .setZones(runtimeAttributes.zones.asJava)
      .setDisks(runtimeAttributes.disks.asJava)
  }
}

object NonPreemptibleJesRuntimeInfo extends JesRuntimeInfoBuilder {
  def apply(commandLine: String, runtimeAttributes: CromwellRuntimeAttributes): JesRuntimeInfo = {
    /*
     It should be impossible for docker to be None here. Enforcing that w/ ADTs seemed more trouble than
      it was worth. If you're ever debugging a NoSuchElementException which leads you here, that means
      the more trouble than worth calculation was incorrect and we should have separate RuntimeAttributes for
      docker and no docker cases
    */
    val dockerImage = runtimeAttributes.docker.get
    val resources = buildResources(runtimeAttributes).setPreemptible(false)
    new NonPreemptibleJesRuntimeInfo(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

object PreemptibleJesRuntimeInfo extends JesRuntimeInfoBuilder {
  def apply(commandLine: String, runtimeAttributes: CromwellRuntimeAttributes): JesRuntimeInfo = {
    // See comment above
    val dockerImage = runtimeAttributes.docker.get
    val resources = buildResources(runtimeAttributes).setPreemptible(true)
    new PreemptibleJesRuntimeInfo(resources, buildDockerExecutor(commandLine, dockerImage))
  }
}

case class NonPreemptibleJesRuntimeInfo(resources: Resources, docker: DockerExecutor) extends JesRuntimeInfo
case class PreemptibleJesRuntimeInfo(resources: Resources, docker: DockerExecutor) extends JesRuntimeInfo