package wom.callable

import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive

/**
  * Parameter documentation quoted from CWL Spec.
  *
  * @param outputPath runtime.outdir: an absolute path to the designated output directory
  * @param tempPath runtime.tmpdir: an absolute path to the designated temporary directory
  * @param cores runtime.cores: number of CPU cores reserved for the tool process
  * @param ram runtime.ram: amount of RAM in mebibytes (2**20) reserved for the tool process
  * @param outputPathSize runtime.outdirSize: reserved storage space available in the designated output directory
  * @param tempPathSize runtime.tmpdirSize: reserved storage space available in the designated temporary directory
  */
case class RuntimeEnvironment(outputPath: String,
                              tempPath: String,
                              cores: Int Refined Positive,
                              ram: Double,
                              outputPathSize: Long,
                              tempPathSize: Long)

