package cromwell.backend.google.pipelines.v2alpha1

import common.util.IntUtil._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineV
import mouse.all._
import org.slf4j.Logger
import squants.information._

object MachineConstraints {
  implicit class EnhancedInformation(val information: Information) extends AnyVal {
    def asMultipleOf(factor: Information): Information = factor * (information / factor).ceil
    def toMBString = information.toString(Megabytes)
    def toMiBString = information.toString(Mebibytes)
  }

  // https://cloud.google.com/compute/docs/instances/creating-instance-with-custom-machine-type
  // https://cloud.google.com/compute/docs/instances/creating-instance-with-custom-machine-type#specifications
  private val minMemoryPerCpu = Gibibytes(0.9)
  private val maxMemoryPerCpu = Gibibytes(6.5)
  private val memoryFactor = Mebibytes(256)

  private def validateCpu(cpu: Int Refined Positive) = cpu.value match {
    // One CPU is cool
    case 1 => 1
    // odd is not allowed, round it up to the next even number
    case odd if odd.isOdd => odd + 1
    case even => even
  }

  private def validateMemory(memory: Information) = memory.asMultipleOf(memoryFactor)

  // Assumes memory and cpu have been validated individually
  private def balanceMemoryAndCpu(memory: Information, cpu: Int) = {
    val memoryPerCpuRatio = memory / cpu.toDouble

    lazy val adjustedMemory = (minMemoryPerCpu * cpu.toDouble) |> validateMemory

    lazy val adjustedCpu = refineV[Positive]((memory / maxMemoryPerCpu).ceil.toInt) match {
      // If for some reason the above yields 0, keep the cpu value unchanged 
      case Left(_) => cpu
      case Right(adjusted) => validateCpu(adjusted)
    }

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
  
  private def logAdjustment(originalCpu: Int, adjustedCpu: Int, originalMemory: Information, adjustedMemory: Information, logger: Logger) = {
    def memoryAdjustmentLog = s"memory was adjusted from ${originalMemory.toMBString} to ${adjustedMemory.toMiBString}"
    def cpuAdjustmentLog = s"cpu was adjusted from $originalCpu to $adjustedCpu"
    
    val message = (originalCpu == adjustedCpu, originalMemory.toMegabytes == adjustedMemory.toMebibytes) match {
      case (true, false) => Option(memoryAdjustmentLog)
      case (false, true) => Option(cpuAdjustmentLog)
      case (false, false) => Option(memoryAdjustmentLog + " and " + cpuAdjustmentLog)
      case _ => None
    }

    message foreach { m => logger.info("To comply with GCE custom machine requirements, " + m) }
  }

  def machineType(memory: Information, cpu: Int Refined Positive, jobLogger: Logger) = {
    val (validCpu, validMemory) = balanceMemoryAndCpu(memory |> validateMemory, cpu |> validateCpu)
    logAdjustment(cpu.value, validCpu, memory, validMemory, jobLogger)
    s"custom-$validCpu-${validMemory.toMebibytes.intValue()}"
  }
}
