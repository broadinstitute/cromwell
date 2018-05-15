package cromwell.backend.google.pipelines.v2alpha1

import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import mouse.all._
import squants.information.{Gigabytes, Information, Megabytes}

object MachineConstraints {
  implicit class EnhancedInformation(val information: Information) extends AnyVal {
    def asMultipleOf(factor: Information): Information = factor * (information / factor).ceil
  }
  
  // https://cloud.google.com/compute/docs/instances/creating-instance-with-custom-machine-type
  // https://cloud.google.com/compute/docs/instances/creating-instance-with-custom-machine-type#specifications
  private val minMemoryPerCpu = Gigabytes(0.9)
  private val maxMemoryPerCpu = Gigabytes(6.5)
  private val memoryFactor = Megabytes(256)

  private def validateCpu(cpu: Int) = cpu match {
    case negative if negative < 1 => 1
    case 1 => 1
    // odd is not allowed, round it up to the next even number
    case odd if odd % 2 == 1 => odd + 1
    // even is fine
    case even => even
  }
  
  private def validateMemory(memory: Information) = memory.asMultipleOf(memoryFactor)

  // Assumes memory and cpu have been validated individually
  private def balanceMemoryAndCpu(memory: Information, cpu: Int) = {
    val memoryPerCpuRatio = memory / cpu.toDouble
    
    lazy val adjustedMemory = (minMemoryPerCpu * cpu.toDouble) |> validateMemory
    // Unfortunately we can't create refined value from non literals, that's why validateCpu takes an Int and not an Int Refined Positive
    lazy val adjustedCpu = (memory / maxMemoryPerCpu).ceil.toInt |> validateCpu

    // If we're under the ratio, top up the memory. Because validMemory will only increase memory (if needed),
    // there's no risk that the call to validMemory will make the ratio invalid
    if (memoryPerCpuRatio < minMemoryPerCpu) {
      cpu -> adjustedMemory
    } else
    // If we're over the ratio, top up the CPU. Because validCpu will only increase CPU (if needed), there's no risk
    // that the call to validCpu will make the ratio invalid
    if (memoryPerCpuRatio > maxMemoryPerCpu) {
      adjustedCpu -> memory
    } else cpu -> memory
  }
  
  def machineType(memory: Information, cpu: Int Refined Positive) = {
    val (validCpu, validMemory) = balanceMemoryAndCpu(memory |> validateMemory, cpu.value |> validateCpu)
    s"custom-$validCpu-${validMemory.toMegabytes.intValue()}"
  }
}
