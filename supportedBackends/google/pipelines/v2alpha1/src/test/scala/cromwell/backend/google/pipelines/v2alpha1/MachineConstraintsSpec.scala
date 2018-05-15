package cromwell.backend.google.pipelines.v2alpha1

import eu.timepit.refined.numeric.Positive
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}
import squants.information.{Gigabytes, Megabytes}
import eu.timepit.refined.refineMV

class MachineConstraintsSpec extends FlatSpec with Matchers {
  behavior of "MachineConstraints"

  it should "generate valid machine types" in {
    val validTypes = Table(
      ("memory", "cpu", "machineTypeString"),
      // Already ok tuple
      (Megabytes(1024), refineMV[Positive](1), "custom-1-1024"),
      // CPU must be even (except if it's 1)
      (Gigabytes(4), refineMV[Positive](3), "custom-4-4096"),
      // Memory must be a multiple of 256
      (Gigabytes(1), refineMV[Positive](1), "custom-1-1024"),
      // Memory / cpu ratio must be > 0.9GB, increase memory
      (Gigabytes(1), refineMV[Positive](4), "custom-4-3840"),
      // Memory / cpu ratio must be < 6.5GB, increase CPU
      (Gigabytes(13.65), refineMV[Positive](1), "custom-4-13824")
    )

    forAll(validTypes) { (memory, cpu, expected) =>
      MachineConstraints.machineType(memory, cpu) shouldBe expected
    }
  }
}
