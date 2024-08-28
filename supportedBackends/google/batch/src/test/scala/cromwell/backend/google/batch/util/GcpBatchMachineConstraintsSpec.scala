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
      ("memory", "cpu", "cpuPlatformOption", "googleLegacyMachineSelection", "machineTypeString"),
      // Already ok tuple
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), None, false, "custom-1-1024"),
      // CPU must be even (except if it's 1)
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), None, false, "custom-4-4096"),
      // Memory must be a multiple of 256
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), None, false, "custom-1-1024"),
      // Memory / cpu ratio must be > 0.9GB, increase memory
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), None, false, "custom-4-3840"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), None, false, "custom-16-14848"),
      // Memory / cpu ratio must be < 6.5GB, increase CPU
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), None, false, "custom-4-14080"),
      // Memory should be an int
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), None, false, "custom-1-1536"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), None, false, "custom-1-1024"),
      // Increase to a cpu selection not valid for n2 below
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, false, "custom-34-31488"),

      // Same tests as above but with legacy machine type selection (cpu and memory as specified. No 'custom machine
      // requirement' adjustments are expected this time, except float->int)
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), None, true, "predefined-1-1024"),
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), None, true, "predefined-3-4096"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), None, true, "predefined-1-1024"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), None, true, "predefined-4-1024"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), None, true, "predefined-16-14336"),
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), None, true, "predefined-1-13977"),
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), None, true, "predefined-1-1520"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), None, true, "predefined-1-1024"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), None, true, "predefined-33-2048"),

      // Same tests but with cascade lake (n2)
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), n2OptionCascadeLake, false, "n2-custom-2-2048"),
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), n2OptionCascadeLake, false, "n2-custom-4-4096"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), n2OptionCascadeLake, false, "n2-custom-2-2048"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), n2OptionCascadeLake, false, "n2-custom-4-4096"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), n2OptionCascadeLake, false, "n2-custom-16-16384"),
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), n2OptionCascadeLake, false, "n2-custom-2-14080"),
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), n2OptionCascadeLake, false, "n2-custom-2-2048"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), n2OptionCascadeLake, false, "n2-custom-2-2048"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), n2OptionCascadeLake, false, "n2-custom-36-36864"),

      // Same tests, but with ice lake. Should produce same results as cascade lake since they're both n2.
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), n2OptionIceLake, false, "n2-custom-2-2048"),
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), n2OptionIceLake, false, "n2-custom-4-4096"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), n2OptionIceLake, false, "n2-custom-2-2048"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), n2OptionIceLake, false, "n2-custom-4-4096"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), n2OptionIceLake, false, "n2-custom-16-16384"),
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), n2OptionIceLake, false, "n2-custom-2-14080"),
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), n2OptionIceLake, false, "n2-custom-2-2048"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), n2OptionIceLake, false, "n2-custom-2-2048"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), n2OptionIceLake, false, "n2-custom-36-36864"),

      // Same tests but with AMD Rome (n2d) #cpu > 16 are in increments of 16
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), n2dOption, false, "n2d-custom-2-1024"),
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), n2dOption, false, "n2d-custom-4-4096"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), n2dOption, false, "n2d-custom-2-1024"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), n2dOption, false, "n2d-custom-4-2048"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), n2dOption, false, "n2d-custom-16-14336"),
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), n2dOption, false, "n2d-custom-2-14080"),
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), n2dOption, false, "n2d-custom-2-1536"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), n2dOption, false, "n2d-custom-2-1024"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](33), n2dOption, false, "n2d-custom-48-24576"),
      (MemorySize(2, MemoryUnit.GB), refineMV[Positive](81), n2dOption, false, "n2d-custom-96-49152"),
      (MemorySize(256, MemoryUnit.GB), refineMV[Positive](128), n2dOption, false, "n2d-custom-96-262144")
    )

    forAll(validTypes) { (memory, cpu, cpuPlatformOption, googleLegacyMachineSelection, expected) =>
      GcpBatchMachineConstraints.machineType(
        memory = memory,
        cpu = cpu,
        cpuPlatformOption = cpuPlatformOption,
        googleLegacyMachineSelection = googleLegacyMachineSelection,
        jobLogger = mock[JobLogger]
      ) shouldBe expected
    }
  }
}
