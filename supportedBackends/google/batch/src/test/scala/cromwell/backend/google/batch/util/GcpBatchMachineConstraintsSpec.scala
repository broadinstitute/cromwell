package cromwell.backend.google.batch.util

import common.assertion.CromwellTimeoutSpec
import common.mock.MockSugar.mock
import cromwell.backend.google.batch.models.GcpBatchRuntimeAttributes
import cromwell.core.logging.JobLogger
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize

class GcpBatchMachineConstraintsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  behavior of "MachineConstraints"

  private val n2OptionCascadeLake = Option(GcpBatchRuntimeAttributes.CpuPlatformIntelCascadeLakeValue)

  private val n2dOption = Option(GcpBatchRuntimeAttributes.CpuPlatformAMDRomeValue)

  private val n2OptionIceLake = Option(GcpBatchRuntimeAttributes.CpuPlatformIntelIceLakeValue)

  it should "generate valid machine types" in {
    val validTypes = Table(
      ("memory", "cpu", "cpuPlatformOption", "standardMachineTypeOption", "googleLegacyMachineSelection", "machineTypeString"),
      // Already ok tuple
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), None, None, false, "custom-1-1024"),
      // CPU must be even (except if it's 1)
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), None, None, false, "custom-4-4096"),
      // Memory must be a multiple of 256
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), None, None, false, "custom-1-1024"),
      // Memory / cpu ratio must be > 0.9GB, increase memory
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), None, None, false, "custom-4-3840"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), None, None, false, "custom-16-14848"),
      // Memory / cpu ratio must be < 6.5GB, increase CPU
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), None, None, false, "custom-4-14080"),
      // Memory should be an int
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), None, None, false, "custom-1-1536"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), None, None, false, "custom-1-1024"),
      // Increase to a cpu selection not valid for n2 below
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, None, false, "custom-34-31488"),

      // Same tests as above but with legacy machine type selection (cpu and memory as specified. No 'custom machine
      // requirement' adjustments are expected this time, except float->int)
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), None, None, true, "predefined-1-1024"),
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), None, None, true, "predefined-3-4096"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), None, None, true, "predefined-1-1024"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), None, None, true, "predefined-4-1024"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), None, None, true, "predefined-16-14336"),
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), None, None, true, "predefined-1-13977"),
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), None, None, true, "predefined-1-1520"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), None, None, true, "predefined-1-1024"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, None, true, "predefined-33-2048"),

      // Same tests but with cascade lake (n2)
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), n2OptionCascadeLake, None, false, "n2-custom-2-2048"),
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), n2OptionCascadeLake, None, false, "n2-custom-4-4096"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), n2OptionCascadeLake, None, false, "n2-custom-2-2048"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), n2OptionCascadeLake, None, false, "n2-custom-4-4096"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), n2OptionCascadeLake, None, false, "n2-custom-16-16384"),
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), n2OptionCascadeLake, None, false, "n2-custom-2-14080"),
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), n2OptionCascadeLake, None, false, "n2-custom-2-2048"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), n2OptionCascadeLake, None, false, "n2-custom-2-2048"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), n2OptionCascadeLake, None, false, "n2-custom-36-36864"),

      // Same tests, but with ice lake. Should produce same results as cascade lake since they're both n2.
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), n2OptionIceLake, None, false, "n2-custom-2-2048"),
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), n2OptionIceLake, None, false, "n2-custom-4-4096"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), n2OptionIceLake, None, false, "n2-custom-2-2048"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), n2OptionIceLake, None, false, "n2-custom-4-4096"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), n2OptionIceLake, None, false, "n2-custom-16-16384"),
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), n2OptionIceLake, None, false, "n2-custom-2-14080"),
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), n2OptionIceLake, None, false, "n2-custom-2-2048"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), n2OptionIceLake, None, false, "n2-custom-2-2048"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), n2OptionIceLake, None, false, "n2-custom-36-36864"),

      // Same tests but with AMD Rome (n2d) #cpu > 16 are in increments of 16
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), n2dOption, None, false, "n2d-custom-2-1024"),
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), n2dOption, None, false, "n2d-custom-4-4096"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), n2dOption, None, false, "n2d-custom-2-1024"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), n2dOption, None, false, "n2d-custom-4-2048"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), n2dOption, None, false, "n2d-custom-16-14336"),
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), n2dOption, None, false, "n2d-custom-2-14080"),
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), n2dOption, None, false, "n2d-custom-2-1536"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), n2dOption, None, false, "n2d-custom-2-1024"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), n2dOption, None, false, "n2d-custom-48-24576"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](81), n2dOption, None, false, "n2d-custom-96-49152"),
      (MemorySize(256, MemoryUnit.GB), refineMV[Positive](128), n2dOption, None, false, "n2d-custom-96-262144"),

      // Test Standard Machine types
      // General-purpose machine family
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("n1-standard-2"), false, "n1-standard-2"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("n1-highmem-2"), false, "n1-highmem-2"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("n1-highcpu-4"), false, "n1-highcpu-4"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("f1-micro"), false, "f1-micro"),

      // Accelerator-optimized machine family
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("a2-highgpu-1g"), false, "a2-highgpu-1g"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("a3-megagpu-8g"), false, "a3-megagpu-8g"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("g2-standard-4"), false, "g2-standard-4"),

      // Other machine families
      // Storage-optimized
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("z3-highmem-88"), false, "z3-highmem-88"),
      // Compute-optimized
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("h3-standard-88"), false, "h3-standard-88"),
      // Memory-optimized
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("m3-ultramem-128"), false, "m3-ultramem-128"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("a2-highgpu-1g"), false, "a2-highgpu-1g"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("a2-highgpu-1g"), false, "a2-highgpu-1g"),

      // Standard machine type overrides legacy selection
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("a2-highgpu-1g"), true, "a2-highgpu-1g"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, Option("a2-highgpu-1g"), false, "a2-highgpu-1g")
    )

    forAll(validTypes) { (memory, cpu, cpuPlatformOption, standardMachineTypeOption, googleLegacyMachineSelection, expected) =>
      GcpBatchMachineConstraints.machineType(
        memory = memory,
        cpu = cpu,
        cpuPlatformOption = cpuPlatformOption,
        standardMachineTypeOption = standardMachineTypeOption,
        googleLegacyMachineSelection = googleLegacyMachineSelection,
        jobLogger = mock[JobLogger]
      ) shouldBe expected
    }
  }
}
