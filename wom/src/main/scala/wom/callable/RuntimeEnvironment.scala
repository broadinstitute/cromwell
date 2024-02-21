package wom.callable

/**
  * Parameter documentation quoted from CWL Spec.
  *
  * @param outputPath runtime.outdir: an absolute path to the designated output directory
  * @param tempPath runtime.tmpdir: an absolute path to the designated temporary directory
  */
case class RuntimeEnvironment(outputPath: String,
                              tempPath: String,
)
