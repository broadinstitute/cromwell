package wom.callable

case class RuntimeEnvironment(outputPath: String,
                              tempPath: String,
                              cores: Int,
                              ram: Double,
                              outputPathSize: Long,
                              tempPathSize: Long)

