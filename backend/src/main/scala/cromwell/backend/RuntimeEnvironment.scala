package cromwell.backend

import java.util.UUID

import cromwell.backend.io.JobPaths
import cromwell.backend.validation.{CpuValidation, MemoryValidation}
import cromwell.core.path.Path
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import wdl4s.parser.MemoryUnit
import wom.callable.RuntimeEnvironment
import wom.format.MemorySize
import wom.values.WomValue

object RuntimeEnvironmentBuilder {

  def apply(runtimeAttributes: Map[String, WomValue], callRoot: Path, callExecutionRoot: Path): MinimumRuntimeSettings => RuntimeEnvironment = {
    minimums =>

      val outputPath: String = callExecutionRoot.pathAsString

      val tempPath: String = {
        val uuid = UUID.randomUUID().toString
        val hash = uuid.substring(0, uuid.indexOf('-'))
        callRoot.resolve(s"tmp.$hash").pathAsString
      }

      val cores: Int Refined Positive = CpuValidation.instanceMin.validate(runtimeAttributes).getOrElse(minimums.cores)

      val memoryInMB: Double =
        MemoryValidation.instance().
          validate(runtimeAttributes).
          map(_.to(MemoryUnit.MB).amount).
          getOrElse(minimums.ram.amount)

      //TODO: Read these from somewhere else
      val outputPathSize: Long = minimums.outputPathSize

      val tempPathSize: Long = minimums.outputPathSize

      RuntimeEnvironment(outputPath, tempPath, cores, memoryInMB, outputPathSize, tempPathSize)
  }

  /**
    * Per the spec:
    *
    * "For cores, ram, outdirSize and tmpdirSize, if an implementation can't provide the actual number of reserved cores
    * during the expression evaluation time, it should report back the minimal requested amount."
    */
  def apply(runtimeAttributes: Map[String, WomValue], jobPaths: JobPaths): MinimumRuntimeSettings => RuntimeEnvironment = {
    this.apply(runtimeAttributes, jobPaths.callRoot, jobPaths.callExecutionRoot)
  }
}

case class MinimumRuntimeSettings(cores: Int Refined Positive = refineMV(1),
                                  ram: MemorySize = MemorySize(4, MemoryUnit.GB),
                                  outputPathSize: Long = Long.MaxValue,
                                  tempPathSize: Long = Long.MaxValue)
