package cromwell.backend.google.pipelines.common

import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import org.slf4j.Logger
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize

object MachineConstraints {
  def machineType(memory: MemorySize,
                  cpu: Int Refined Positive,
                  cpuPlatformOption: Option[String],
                  googleLegacyMachineSelection: Boolean,
                  jobLogger: Logger
  ): String =
    if (googleLegacyMachineSelection) {
      s"predefined-$cpu-${memory.to(MemoryUnit.MB).amount.intValue()}"
    } else {
      // Users specify a CPU platform in their WDL, but GCP also needs to know which machine type to use.
      // The below logic infers the machine type from the requested CPU.
      // The heuristic we're using is: find the newest 'General Purpose' type that supports the given CPU.
      // https://cloud.google.com/compute/docs/machine-resource
      // For example, if someone requests Intel Cascade Lake as their CPU platform, then infer the n2 machine type.
      // Infer n2d from AMD Rome, etc.
      val customMachineType =
        cpuPlatformOption match {
          case Some(PipelinesApiRuntimeAttributes.CpuPlatformIntelIceLakeValue) => N2CustomMachineType
          case Some(PipelinesApiRuntimeAttributes.CpuPlatformIntelCascadeLakeValue) => N2CustomMachineType
          case Some(PipelinesApiRuntimeAttributes.CpuPlatformAMDRomeValue) => N2DCustomMachineType
          case _ => N1CustomMachineType
        }
      customMachineType.machineType(memory, cpu, jobLogger)
    }
}
