package cromwell.binding

import com.google.api.services.genomics.model.Disk
import cromwell.binding.AstTools.{AstNodeName, EnhancedAstNode}
import cromwell.binding.RuntimeAttributes.{Defaults, LocalDisk}
import cromwell.parser.RuntimeKey._
import cromwell.parser.WdlParser.{Ast, AstList}
import cromwell.parser.{BackendType, MemorySize, RuntimeKey}
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.language.postfixOps

object RuntimeAttributes {
  /** Convenience class to convert to Google `Disk` class. */
  private case class LocalDisk(name: String, sizeGb: Long, diskType: String) {
    def toDisk: Disk = new Disk().setSizeGb(sizeGb).setType(diskType).setName(name)
  }

  private val log = LoggerFactory.getLogger("RuntimeAttributes")

  /** Fallback values if values for these keys are not specified in a task "runtime" stanza. */
  object Defaults {
    val Cpu = 1
    val Disk = Seq(LocalDisk("local-disk", 100L, "LOCAL_SSD").toDisk)
    val FailOnStderr = false
    val MemoryInBytes = MemorySize.GB.toBytes(2)
    val Preemptible = false
    val Zones = Array("us-central1-a")
  }

  def apply(ast: Ast, backendType: BackendType): RuntimeAttributes = {
    def processRuntimeAttribute(ast: Ast): (String, String) = {
      (ast.getAttribute("key").sourceString, ast.getAttribute("value").sourceString)
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
    val actualKeys = attributeMap.keySet
    val knownKeys = RuntimeKey.values().map { _.key }

    // Finding keys unknown to any backend is an error.
    val unknownKeys = actualKeys -- knownKeys
    if (unknownKeys.nonEmpty) {
      throw new RuntimeException(s"Unknown keys found in runtime configuration: "
        + unknownKeys.toSeq.sorted.mkString(", "))
    }

    // Missing required keys for this backend is an error.
    val requiredKeys = for {
      key <- RuntimeKey.values().toSet
      if key.isMandatory(backendType)
    } yield key.key

    val missingRequiredKeys = requiredKeys -- actualKeys
    if (missingRequiredKeys.nonEmpty) {
      throw new RuntimeException(s"Missing required keys in runtime configuration for backend '$backendType': "
        + missingRequiredKeys.toSeq.sorted.mkString(", "))
    }

    // Keys supported on this backend.
    val supportedKeys = for {
      key <- RuntimeKey.values()
      if key.supports(backendType)
    } yield key.key

    // Warn if keys are found that are unsupported on this backend.  This is not necessarily an error, at this point
    // the keys are known to be supported on at least one backend other than this.
    val unsupportedKeys = actualKeys -- supportedKeys
    val warnings = if (unsupportedKeys.isEmpty) {
      Seq.empty
    } else {
      Seq(s"Found unsupported keys for backend '$backendType': " + unsupportedKeys.toSeq.sorted.mkString(", "))
    }
    warnings foreach log.warn
    RuntimeAttributes(attributeMap, warnings)
  }
}

/**
 * Parses configuration values per WDL/JES specs.  Constructor is private to force use of companion object
 * `apply` method.
 */
case class RuntimeAttributes private(attributes: Map[String, String], warnings: Seq[String]) {

  private def attr(key: RuntimeKey): Option[String] = attributes.get(key.key)

  val docker: Option[String] = attr(DOCKER)
  val failOnStderr: Boolean = attr(FAIL_ON_STDERR) map { _.toBoolean } getOrElse Defaults.FailOnStderr
  val cpu: Long = attr(CPU) map { _.toLong } getOrElse Defaults.Cpu
  val preemptible: Boolean = attr(PREEMPTIBLE) map { _.toBoolean } getOrElse Defaults.Preemptible
  val defaultZones: Seq[String] = (attr(DEFAULT_ZONES) map { _.split("\\s+")} getOrElse Defaults.Zones).toSeq

  val defaultDisks: Seq[Disk] = {
    val taskSpecifiedAttributes =
      for {
      // This toSeq on the Option creates a Seq with zero or one element.  If there are zero elements,
      // nothing will be returned from the overall for comprehension.
        fullString <- attr(DEFAULT_DISKS).toSeq
        diskString <- fullString.split(",\\s*")
        Array(name, sizeGb, diskType) = diskString.split("\\s+")
      } yield LocalDisk(name, sizeGb.toLong, diskType).toDisk

    if (taskSpecifiedAttributes.isEmpty) Defaults.Disk else taskSpecifiedAttributes
  }

  val memoryGB: Double = {
    val taskSpecifiedBytes = {
      for {
        memory <- attr(MEMORY).toSeq
        Array(prefix, suffix) = memory.split("\\s+")
        memorySize <- MemorySize.values().find(_.suffix == suffix)
      } yield memorySize.toBytes(prefix.toDouble)
    }
    val memoryInBytes = taskSpecifiedBytes.headOption.getOrElse(Defaults.MemoryInBytes)
    MemorySize.GB.fromBytes(memoryInBytes)
  }
}
