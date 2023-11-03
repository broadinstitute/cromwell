package cromwell.backend.google.pipelines.common

import cromwell.backend.google.pipelines.common.CustomMachineType._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineV
import mouse.all._
import org.slf4j.Logger
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import math.{log, pow}

case class PredefinedMachineType(cpuCount: Int, gcpString: String)
object PredefinedMachineType {
  val c3Standard_4 = PredefinedMachineType(4, "c3-standard-4")
  val c3Standard_8 = PredefinedMachineType(8, "c3-standard-8")
  val c3Standard_22 = PredefinedMachineType(22, "c3-standard-22")
  val c3Standard_44 = PredefinedMachineType(44, "c3-standard-44")
  val c3Standard_88 = PredefinedMachineType(88, "c3-standard-88")
  val c3Standard_176 = PredefinedMachineType(176, "c3-standard-176")

  def getClosestC3Machine(requestedMemory: MemorySize, requestedCpu: Refined[Int, Positive],  jobLogger: Logger): String = {
    val adjustedMemory: MemorySize = MemorySize(requestedCpu.value * 4.0, MemoryUnit.GB)
    if (adjustedMemory != requestedMemory) {
      jobLogger.info(s"Adjusting memory from ${requestedMemory.amount} to ${adjustedMemory.amount} in order to match GCP requirements for the requested CPU.")
    }
    val machine = requestedCpu.value match {
      case cpu if cpu <= 4 => c3Standard_4
      case cpu if cpu > 4 && cpu <= 8 => c3Standard_8
      case cpu if cpu > 8 && cpu <= 22 => c3Standard_22
      case cpu if cpu > 22 && cpu <= 44 => c3Standard_44
      case cpu if cpu > 44 && cpu <= 88 => c3Standard_88
      case cpu if cpu > 88 => c3Standard_176
    }

    if(machine.cpuCount != requestedCpu.value) {
      jobLogger.info(s"Rounding up CPU count from ${requestedCpu} to ${machine.cpuCount} in order to match GCP requirements for the requested CPU.")
    }
    machine.gcpString
  }
}
