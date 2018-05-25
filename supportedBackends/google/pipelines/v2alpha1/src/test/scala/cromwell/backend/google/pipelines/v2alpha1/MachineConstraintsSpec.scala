package cromwell.backend.google.pipelines.v2alpha1

import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}
import org.slf4j.helpers.NOPLogger
import squants.information.{Gigabytes, Megabytes}

class MachineConstraintsSpec extends FlatSpec with Matchers {
  behavior of "MachineConstraints"

  it should "generate valid machine types" in {
    val validTypes = Table(
      ("memory", "cpu", "machineTypeString"),
      // Already ok tuple
      (Megabytes(1024), refineMV[Positive](1), "custom-1-1024"),
      // CPU must be even (except if it's 1)
      (Gigabytes(4), refineMV[Positive](3), "custom-4-3840"),
      // Memory must be a multiple of 256
      (Gigabytes(1), refineMV[Positive](1), "custom-1-1024"),
      // Memory / cpu ratio must be > 0.9GB, increase memory
      (Gigabytes(1), refineMV[Positive](4), "custom-4-3840"),
      (Gigabytes(14), refineMV[Positive](16), "custom-16-14848"),
      // Memory / cpu ratio must be < 6.5GB, increase CPU
      (Gigabytes(13.65), refineMV[Positive](1), "custom-2-13056"),
      // Memory should be an int
      (Megabytes(1520.96), refineMV[Positive](1), "custom-1-1536"),
      (Megabytes(1024.0), refineMV[Positive](1), "custom-1-1024")
    )

    forAll(validTypes) { (memory, cpu, expected) =>
      MachineConstraints.machineType(memory, cpu, NOPLogger.NOP_LOGGER) shouldBe expected
    }
  }
}
