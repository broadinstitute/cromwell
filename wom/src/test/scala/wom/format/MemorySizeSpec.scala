package wom.format

import wdl4s.parser.MemoryUnit
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class MemorySizeSpec extends AnyFlatSpec with Matchers {
  behavior of "MemorySize"
  
  it should "provide the rounded up multiple of a number" in {
    MemorySize(2, MemoryUnit.GB).to(MemoryUnit.MB).asRoundedUpMultipleOf(256) shouldBe MemorySize(2048, MemoryUnit.MB)
    MemorySize(2000, MemoryUnit.MB).asRoundedUpMultipleOf(256) shouldBe MemorySize(2048, MemoryUnit.MB)
    MemorySize(3.75, MemoryUnit.GB).to(MemoryUnit.MB).asRoundedUpMultipleOf(256) shouldBe MemorySize(3840, MemoryUnit.MB)
    assertThrows[IllegalArgumentException](MemorySize(2000, MemoryUnit.MB).asRoundedUpMultipleOf(0))
    assertThrows[IllegalArgumentException](MemorySize(2000, MemoryUnit.MB).asRoundedUpMultipleOf(-5))
  }
}
