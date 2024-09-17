package cromwell.backend.google.batch.util

import cromwell.backend.google.batch.models.{
  GcpBatchRuntimeAttributes,
  N1CustomMachineType,
  N2CustomMachineType,
  N2DCustomMachineType,
  StandardMachineType
}
import cromwell.core.logging.JobLogger
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import scala.util.matching.Regex
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize

object GcpBatchMachineConstraints {
  private val machineTypePattern: Regex = """^\w{2}\w?-\w+-\w+$""".r

  def machineType(memory: MemorySize,
                  cpu: Int Refined Positive,
                  cpuPlatformOption: Option[String],
                  googleLegacyMachineSelection: Boolean,
                  jobLogger: JobLogger
  ): String =
    if (isStandardMachineType(cpuPlatformOption.getOrElse(""))) {
      StandardMachineType(cpuPlatformOption.getOrElse("")).machineType
    } else if (googleLegacyMachineSelection) {
      s"predefined-$cpu-${memory.to(MemoryUnit.MB).amount.intValue()}"
    } else {
      // If someone requests Intel Cascade Lake or Intel Ice Lake as their CPU platform then switch the machine type to n2.
      // Similarly, CPU platform of AMD Rome corresponds to the machine type n2d.
      val customMachineType =
        cpuPlatformOption match {
          case Some(GcpBatchRuntimeAttributes.CpuPlatformIntelCascadeLakeValue) => N2CustomMachineType
          case Some(GcpBatchRuntimeAttributes.CpuPlatformIntelIceLakeValue) => N2CustomMachineType
          case Some(GcpBatchRuntimeAttributes.CpuPlatformAMDRomeValue) => N2DCustomMachineType
          case _ => N1CustomMachineType
        }
      customMachineType.machineType(memory, cpu, jobLogger)
    }

  def isStandardMachineType(machineType: String): Boolean = machineTypePattern.matches(machineType)
}
