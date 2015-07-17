package cromwell.binding

import com.google.api.services.genomics.model.Disk
import cromwell.binding.AstTools.{AstNodeName, EnhancedAstNode}
import cromwell.binding.RuntimeAttributes.{Defaults, LocalDisk}
import cromwell.parser.MemorySize
import cromwell.parser.WdlParser.{Ast, AstList}

import scala.collection.JavaConverters._
import scala.language.postfixOps

object RuntimeAttributes {
  /** Convenience class to convert to Google `Disk` class. */
  private case class LocalDisk(name: String, sizeGb: Long, diskType: String) {
    def toDisk: Disk = new Disk().setSizeGb(sizeGb).setType(diskType).setName(name)
  }

  /** Fallback values if values for these keys are not specified in a task "runtime" stanza. */
  object Defaults {
    val Cpu = 1
    val Disk = Seq(LocalDisk("local-disk", 100L, "LOCAL_SSD").toDisk)
    val FailOnStderr = false
    val MemoryInBytes = MemorySize.GB.toBytes(2)
    val Preemptible = false
    val Zones = Array("us-central1-a")
  }

  def apply(ast: Ast): RuntimeAttributes = {
    def processRuntimeAttribute(ast: Ast): (String, String) = {
      (ast.getAttribute("key").sourceString(), ast.getAttribute("value").sourceString())
    }

    def processRuntimeAttributes(astList: AstList): Map[String, String] = {
      astList.asScala.toVector map { a => processRuntimeAttribute(a.asInstanceOf[Ast]) } toMap
    }

    val asts = ast.findAsts(AstNodeName.Runtime)
    if (asts.size > 1) throw new UnsupportedOperationException("Only one runtime block may be defined per task")
    val astList = asts.headOption map {
      _.getAttribute("map").asInstanceOf[AstList]
    }
    val attributeMap = astList map processRuntimeAttributes getOrElse Map.empty[String, String]

    RuntimeAttributes(attributeMap)
  }
}

/**
 * Parses configuration values per WDL/JES specs.  Constructor is private to force use of companion object
 * `apply` method.
 */
case class RuntimeAttributes private(attributes: Map[String, String]) {

  def attribute(attr: String): Option[String] = attributes.get(attr)

  def docker: Option[String] = attribute("docker")
  def failOnStderr: Boolean = attribute("failOnStderr") map { _.toBoolean } getOrElse Defaults.FailOnStderr
  def cpu: Long = attribute("cpu") map { _.toLong } getOrElse Defaults.Cpu
  def preemptible: Boolean = attribute("preemptible") map { _.toBoolean } getOrElse Defaults.Preemptible
  def defaultZones: Seq[String] = (attribute("defaultZones") map { _.split("\\s+")} getOrElse Defaults.Zones).toSeq

  def defaultDisks: Seq[Disk] = {
    val taskSpecifiedAttributes =
      for {
      // This toSeq on the Option creates a Seq with zero or one element.  If there are zero elements,
      // nothing will be returned from the overall for comprehension.
        fullString <- attribute("defaultDisks").toSeq
        diskString <- fullString.split(",\\s*")
        Array(name, sizeGb, diskType) = diskString.split("\\s+")
      } yield LocalDisk(name, sizeGb.toLong, diskType).toDisk

    if (taskSpecifiedAttributes.isEmpty) Defaults.Disk else taskSpecifiedAttributes
  }

  def memoryGB: Double = {
    val taskSpecifiedBytes = {
      for {
        memory <- attribute("memory").toSeq
        Array(prefix, suffix) = memory.split("\\s+")
        memorySize <- MemorySize.values().find(_.suffix == suffix)
      } yield memorySize.toBytes(prefix.toDouble)
    }
    val memoryInBytes = taskSpecifiedBytes.headOption.getOrElse(Defaults.MemoryInBytes)
    MemorySize.GB.fromBytes(memoryInBytes)
  }
}
