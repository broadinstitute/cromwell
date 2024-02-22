package cromwell.backend

import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize

case class MinimumRuntimeSettings(cores: Int Refined Positive = refineMV(1),
                                  ram: MemorySize = MemorySize(4, MemoryUnit.GB),
                                  outputPathSize: Long = Long.MaxValue,
                                  tempPathSize: Long = Long.MaxValue
)
