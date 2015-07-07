package cromwell.engine.backend.jes

import scala.collection.JavaConverters._
import com.google.api.services.genomics.model.{Disk, DockerExecutor, Resources}
import cromwell.binding.Call

object JesRuntimeInfo {
  val DefaultCpu = 1L
  val DefaultMemory = 2.toDouble
  val DefaultZones = Vector("us-central1-a")
  val DefaultPreemptible = false
  val DefaultDisks = Vector(LocalDisk("local-disk", 100L, "LOCAL_SSD"))

  def apply(commandLine: String, call: Call): JesRuntimeInfo = {
    val dockerImage = call.docker.getOrElse(throw new IllegalArgumentException("FIXME: WHAT TO DO?"))

    new JesRuntimeInfo(buildResources(call), buildDockerExecutor(commandLine, dockerImage))
  }

  def buildDockerExecutor(commandLine: String, dockerImage: String): DockerExecutor = {
    val docker = new DockerExecutor()
    docker.setImage(dockerImage).setCmd(s"/bin/bash -c '$commandLine'")
  }

  def buildResources(call: Call): Resources = {
    new Resources()
      .setRamGb(DefaultMemory) // FIXME: We'll need to parse the task.runtimeAttributes (probably in binding/package.scala) to convert to GB
      .setCpu(DefaultCpu)
      .setPreemptible(DefaultPreemptible)
      .setZones(DefaultZones.asJava)
      .setDisks(DefaultDisks.map{_.toDisk}.asJava)
  }

  case class LocalDisk(name: String, sizeGb: Long, diskType: String) {
    def toDisk: Disk = {
      val localDisk = new Disk()
      localDisk.setSizeGb(sizeGb).setType(diskType).setName(name)
    }
  }
}

case class JesRuntimeInfo(resources: Resources, docker: DockerExecutor)
