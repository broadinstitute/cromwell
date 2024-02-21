package cromwell.backend

import java.util.UUID

import cromwell.backend.io.JobPaths
import cromwell.core.path.Path
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import wdl4s.parser.MemoryUnit
import wom.callable.RuntimeEnvironment
import wom.format.MemorySize

object RuntimeEnvironmentBuilder {

  def apply(callRoot: Path,
            callExecutionRoot: Path
  ): MinimumRuntimeSettings => RuntimeEnvironment = { minimums =>
    val outputPath: String = callExecutionRoot.pathAsString

    val tempPath: String = {
      val uuid = UUID.randomUUID().toString
      val hash = uuid.substring(0, uuid.indexOf('-'))
      callRoot.resolve(s"tmp.$hash").pathAsString
    }

    RuntimeEnvironment(outputPath, tempPath)
  }

  /**
    * Per the spec:
    *
    * "For cores, ram, outdirSize and tmpdirSize, if an implementation can't provide the actual number of reserved cores
    * during the expression evaluation time, it should report back the minimal requested amount."
    */
  def apply(jobPaths: JobPaths): MinimumRuntimeSettings => RuntimeEnvironment =
    this.apply(jobPaths.callRoot, jobPaths.callExecutionRoot)
}

case class MinimumRuntimeSettings(cores: Int Refined Positive = refineMV(1),
                                  ram: MemorySize = MemorySize(4, MemoryUnit.GB),
                                  outputPathSize: Long = Long.MaxValue,
                                  tempPathSize: Long = Long.MaxValue
)
