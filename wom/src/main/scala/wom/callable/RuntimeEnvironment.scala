package wom.callable

/**
  * from CWL Spec:
  * runtime.outdir: an absolute path to the designated output directory
  * runtime.tmpdir: an absolute path to the designated temporary directory
  * runtime.cores: number of CPU cores reserved for the tool process
  * runtime.ram: amount of RAM in mebibytes (2**20) reserved for the tool process
  * runtime.outdirSize: reserved storage space available in the designated output directory
  * runtime.tmpdirSize: reserved storage space available in the designated temporary directory
  */
case class RuntimeEnvironment(outputPath: String,
                              tempPath: String,
                              cores: Int,
                              ram: Int,
                              outputPathSize: Int,
                              tempPathSize: Int
                             )

