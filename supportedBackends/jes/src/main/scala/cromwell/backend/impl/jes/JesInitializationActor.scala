package cromwell.backend.impl.jes

import akka.actor.Props
import cromwell.backend.BackendLifecycleActor.WorkflowAbortResponse
import cromwell.backend.impl.jes.JesInitializationActor._
import cromwell.backend.impl.jes.io.{JesAttachedDisk, JesWorkingDisk}
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core._
import wdl4s.Call
import wdl4s.types.{WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlArray, WdlInteger, WdlString, WdlValue}

import scala.concurrent.Future
import scalaz.Scalaz._

object JesInitializationActor {
  //TODO: check if these need to be configurable.
  val CpuDefaultValue = 1
  val ZonesKey = "zones"
  val ZoneDefaultValue = Seq("us-central1-a")
  val PreemptibleKey = "preemptible"
  val PreemptibleDefaultValue = 0
  val BootDiskSizeKey = "bootDiskSizeGb"
  val BootDiskSizeDefaultValue = 10
  val MemoryDefaultValue = "2 GB"
  val DisksKey = "disks"
  val DisksDefaultValue = s"${JesWorkingDisk.Name} 10 SSD"
  val DefaultJesWorkingDisk = JesAttachedDisk.parse(DisksDefaultValue).get
  val FailOnStderrDefaultValue = false
  val ContinueOnRcDefaultValue = 0

  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], configurationDescriptor: BackendConfigurationDescriptor): Props =
    Props(new JesInitializationActor(workflowDescriptor, calls, configurationDescriptor))
}

class JesInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                             override val calls: Seq[Call],
                             override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendWorkflowInitializationActor {
  /**
    * Abort all initializations.
    */
  override def abortInitialization(): Future[WorkflowAbortResponse] = ???

  //TODO: Workflow options may need to be validated for JES.

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Unit] = Future.successful(())

  /**
    * Validates runtime attributes for one specific call.
    *
    * @param runtimeAttributes Runtime Attributes with already evaluated values.
    * @return If all entries from runtime attributes section are valid Success otherwise
    *         Failure with the aggregation of errors.
    */
  override def validateRuntimeAttributes(runtimeAttributes: EvaluatedRuntimeAttributes): Future[ErrorOr[Unit]] = {
    Future {
      val cpu = validateCpu(runtimeAttributes.get(Cpu), CpuDefaultValue.successNel)
      val zones = validateZone(runtimeAttributes.get(ZonesKey))
      val preemptible = validatePreemptible(runtimeAttributes.get(PreemptibleKey))
      val bootDiskSize = validateBootDisk(runtimeAttributes.get(BootDiskSizeKey))
      val memory = validateMemory(runtimeAttributes.get(Memory), parseMemoryString(WdlString(MemoryDefaultValue)))
      val disks = validateLocalDisks(runtimeAttributes.get(DisksKey))
      val docker = validateDocker(runtimeAttributes.get(Docker), "Failed to get Docker mandatory key from runtime attributes".failureNel)
      val failOnStderr = validateFailOnStderr(runtimeAttributes.get(FailOnStderr), FailOnStderrDefaultValue.successNel)
      val continueOnReturnCode = validateContinueOnReturnCode(runtimeAttributes.get(ContinueOnReturnCode),
        ContinueOnReturnCodeSet(Set(ContinueOnRcDefaultValue)).successNel)
      (cpu |@| zones |@| preemptible |@| bootDiskSize |@| memory |@| disks |@| docker |@| failOnStderr |@| continueOnReturnCode) {
        (_, _, _, _, _, _, _, _, _)
      }
    }
  }

  private def validateZone(zoneValue: Option[WdlValue]): ErrorOr[Vector[String]] = {
    zoneValue match {
      case Some(WdlString(s)) => s.split("\\s+").toVector.successNel
      case Some(WdlArray(wdlType, value)) if wdlType.memberType == WdlStringType =>
        value.map(_.valueString).toVector.successNel
      case Some(_) => s"Expecting $ZonesKey runtime attribute to be either a whitespace separated String or an Array[String]".failureNel
      case None => ZoneDefaultValue.toVector.successNel
    }
  }

  private def validatePreemptible(preemptible: Option[WdlValue]): ErrorOr[Int] = {
    val preemptibleValidation = preemptible.map(validateInt).getOrElse(PreemptibleDefaultValue.successNel)
    if (preemptibleValidation.isFailure) {
      s"Expecting $PreemptibleKey runtime attribute to be an Integer".failureNel
    }
    else {
      preemptibleValidation
    }
  }

  private def validateBootDisk(diskSize: Option[WdlValue]): ErrorOr[Int] = diskSize match {
    case Some(x) if WdlIntegerType.isCoerceableFrom(x.wdlType) =>
      WdlIntegerType.coerceRawValue(x) match {
        case scala.util.Success(x: WdlInteger) => x.value.intValue.successNel
        case scala.util.Success(unhandled) => s"Coercion was expected to create an Integer but instead got $unhandled".failureNel
        case scala.util.Failure(t) => s"Expecting $BootDiskSizeKey runtime attribute to be an Integer".failureNel
      }
    case None => BootDiskSizeDefaultValue.successNel
  }

  private def validateLocalDisks(value: Option[WdlValue]): ErrorOr[Seq[JesAttachedDisk]] = {
    val nels = value match {
      case Some(WdlString(s)) => s.split(",\\s*").toSeq.map(validateLocalDisk)
      case Some(WdlArray(wdlType, seq)) if wdlType.memberType == WdlStringType =>
        seq.map(_.valueString).map(validateLocalDisk)
      case Some(_) =>
        Seq(s"Expecting $DisksKey runtime attribute to be a comma separated String or Array[String]".failureNel[JesAttachedDisk])
      case None => Seq(DefaultJesWorkingDisk).map(_.successNel)
    }

    val emptyDiskNel = Vector.empty[JesAttachedDisk].successNel[String]
    val disksNel = nels.foldLeft(emptyDiskNel)((acc, v) => (acc |@| v) { (a, v) => a :+ v })

    disksNel map {
      case disks if disks.exists(_.name == JesWorkingDisk.Name) => disks
      case disks => disks :+ DefaultJesWorkingDisk
    }
  }

  private def validateLocalDisk(disk: String): ErrorOr[JesAttachedDisk] = {
    JesAttachedDisk.parse(disk) match {
      case scala.util.Success(localDisk) => localDisk.successNel
      case scala.util.Failure(ex) => ex.getMessage.failureNel
    }
  }
}
