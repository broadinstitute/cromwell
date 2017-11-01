package cromwell.backend

import cromwell.backend.io.JobPaths
import cromwell.backend.validation.{CpuValidation, MemoryValidation}
import wdl4s.parser.MemoryUnit
import wom.callable.RuntimeEnvironment
import wom.values.WomValue

object RuntimeEnvironmentBuilder {

  /**
    * Per the spec:
    *
    * "For cores, ram, outdirSize and tmpdirSize, if an implementation can't provide the actual number of reserved cores
    * during the expression evaluation time, it should report back the minimal requested amount."
    */
   def apply(runtimeAttributes: Map[String, WomValue], jobPaths: JobPaths): MinimumRuntimeEnvironment => RuntimeEnvironment = {
     minimums =>

       val outputPath = jobPaths.callRoot

       val tempPath = jobPaths.callRoot

        val cores: Int = CpuValidation.instance.validate(runtimeAttributes).getOrElse(minimums.cores)

        val memoryInMiB: Double =
          MemoryValidation.instance().
            validate(runtimeAttributes).
            map(_.to(MemoryUnit.MiB).amount).
            getOrElse(minimums.ram.amount)

       //TODO: Read these from somewhere else
       val outputPathSize: Long = minimums.outputPathSize

       val tempPathSize: Long = minimums.outputPathSize

       RuntimeEnvironment(outputPath.pathAsString, tempPath.pathAsString, cores, memoryInMiB, outputPathSize, tempPathSize)
  }
}

case class MinimumRuntimeEnvironment(cores: Int = 1, ram: MemorySize = MemorySize(4, MemoryUnit.GiB), outputPathSize: Long = Long.MaxValue, tempPathSize: Long = Long.MaxValue)
