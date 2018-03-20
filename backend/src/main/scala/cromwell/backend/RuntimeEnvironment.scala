package cromwell.backend

import java.util.UUID

import cromwell.backend.io.JobPaths
import cromwell.backend.validation.{CpuValidation, MemoryValidation}
import wdl4s.parser.MemoryUnit
import wom.callable.RuntimeEnvironment
import wom.format.MemorySize
import wom.values.WomValue

object RuntimeEnvironmentBuilder {

  /**
    * Per the spec:
    *
    * "For cores, ram, outdirSize and tmpdirSize, if an implementation can't provide the actual number of reserved cores
    * during the expression evaluation time, it should report back the minimal requested amount."
    */
   def apply(runtimeAttributes: Map[String, WomValue], jobPaths: JobPaths): MinimumRuntimeSettings => RuntimeEnvironment = {
     minimums =>

       val outputPath: String = jobPaths.callExecutionRoot.pathAsString

       val tempPath: String = {
         val uuid = UUID.randomUUID().toString
         val hash = uuid.substring(0, uuid.indexOf('-'))
         jobPaths.callRoot.resolve(s"tmp.$hash").pathAsString
       }

       val cores: Int = CpuValidation.instanceMin.validate(runtimeAttributes).getOrElse(minimums.cores)

       val memoryInMiB: Double =
         MemoryValidation.instance().
           validate(runtimeAttributes).
           map(_.to(MemoryUnit.MiB).amount).
           getOrElse(minimums.ram.amount)

       //TODO: Read these from somewhere else
       val outputPathSize: Long = minimums.outputPathSize

       val tempPathSize: Long = minimums.outputPathSize

       RuntimeEnvironment(outputPath, tempPath, cores, memoryInMiB, outputPathSize, tempPathSize)
  }
}

case class MinimumRuntimeSettings(cores: Int = 1,
                                  ram: MemorySize = MemorySize(4, MemoryUnit.GiB),
                                  outputPathSize: Long = Long.MaxValue,
                                  tempPathSize: Long = Long.MaxValue)
