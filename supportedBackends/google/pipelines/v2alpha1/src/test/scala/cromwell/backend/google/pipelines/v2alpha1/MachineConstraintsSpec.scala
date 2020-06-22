package cromwell.backend.google.pipelines.v2alpha1

import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}
import org.slf4j.helpers.NOPLogger
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize

class MachineConstraintsSpec extends FlatSpec with Matchers {
  behavior of "MachineConstraints"

  it should "generate valid machine types" in {
    val validTypes = Table(
      ("memory", "cpu", "googleLegacyMachineSelection", "machineTypeString"),
      // Already ok tuple
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), false, "custom-1-1024"),
      // CPU must be even (except if it's 1)
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), false, "custom-4-4096"),
      // Memory must be a multiple of 256
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), false, "custom-1-1024"),
      // Memory / cpu ratio must be > 0.9GB, increase memory
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), false, "custom-4-3840"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), false, "custom-16-14848"),
      // Memory / cpu ratio must be < 6.5GB, increase CPU
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), false, "custom-4-14080"),
      // Memory should be an int
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), false, "custom-1-1536"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), false, "custom-1-1024"),

      // Same tests as above but with legacy machine type selection (cpu and memory as specified. No 'custom machine
      // requirement' adjustments are expected this time, except float->int)

      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), true, "predefined-1-1024"),
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), true, "predefined-3-4096"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), true, "predefined-1-1024"),
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), true, "predefined-4-1024"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), true, "predefined-16-14336"),
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), true, "predefined-1-13977"),
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), true, "predefined-1-1520"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), true, "predefined-1-1024"),
    )

    forAll(validTypes) { (memory, cpu, googleLegacyMachineSelection, expected) =>
      MachineConstraints.machineType(memory, cpu, googleLegacyMachineSelection, NOPLogger.NOP_LOGGER) shouldBe expected
    }
  }
}
