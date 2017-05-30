package cromwell.backend

import com.typesafe.config.ConfigFactory

object TestConfig {

  lazy val sampleBackendRuntimeConfigString =
    s"""
       |default-runtime-attributes {
       |    failOnStderr: false
       |    continueOnReturnCode: 0
       |    memory: "2 GB"
       |}
       |""".stripMargin

  lazy val allBackendRuntimeAttrsString =
    s"""
       |default-runtime-attributes {
       |  cpu: 1
       |  failOnStderr: false
       |  continueOnReturnCode: 0
       |  memory: "1 GB"
       |  bootDiskSizeGb: 10
       |  disks: "local-disk 10 SSD"
       |  noAddress: false
       |  preemptible: 0
       |  zones: ["us-central1-a", "us-central1-b"]
       |}
     """.stripMargin

  lazy val sampleBackendRuntimeConfig = ConfigFactory.parseString(sampleBackendRuntimeConfigString)

  lazy val allRuntimeAttrsConfig = ConfigFactory.parseString(allBackendRuntimeAttrsString).getConfig("default-runtime-attributes")

  lazy val optionalRuntimeConfig = sampleBackendRuntimeConfig.getConfig("default-runtime-attributes")

  lazy val globalConfig = ConfigFactory.load()

  lazy val emptyConfig = ConfigFactory.empty()

  lazy val emptyBackendConfigDescriptor = BackendConfigurationDescriptor(emptyConfig, globalConfig)

  lazy val backendRuntimeConfigDescriptor = BackendConfigurationDescriptor(sampleBackendRuntimeConfig, emptyConfig)
}
