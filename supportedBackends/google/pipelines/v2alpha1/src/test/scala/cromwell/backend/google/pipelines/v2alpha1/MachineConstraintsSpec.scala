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
      ("memory", "cpu", "machineTypeString"),
      // Already ok tuple
      (MemorySize(1024, MemoryUnit.MB), refineMV[Positive](1), "custom-1-1024"),
      // CPU must be even (except if it's 1)
      (MemorySize(4, MemoryUnit.GB), refineMV[Positive](3), "custom-4-4096"),
      // Memory must be a multiple of 256
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](1), "custom-1-1024"),
      // Memory / cpu ratio must be > 0.9GB, increase memory
      (MemorySize(1, MemoryUnit.GB), refineMV[Positive](4), "custom-4-3840"),
      (MemorySize(14, MemoryUnit.GB), refineMV[Positive](16), "custom-16-14848"),
      // Memory / cpu ratio must be < 6.5GB, increase CPU
      (MemorySize(13.65, MemoryUnit.GB), refineMV[Positive](1), "custom-4-14080"),
      // Memory should be an int
      (MemorySize(1520.96, MemoryUnit.MB), refineMV[Positive](1), "custom-1-1536"),
      (MemorySize(1024.0, MemoryUnit.MB), refineMV[Positive](1), "custom-1-1024")
    )

    forAll(validTypes) { (memory, cpu, expected) =>
      MachineConstraints.machineType(memory, cpu, NOPLogger.NOP_LOGGER) shouldBe expected
    }
  }
}
