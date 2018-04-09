package wom.format

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.parser.MemoryUnit

class MemorySizeSpec extends FlatSpec with Matchers {
  behavior of "MemorySize"
  
  it should "provide the rounded up multiple of a number" in {
    MemorySize(2000, MemoryUnit.MB).asRoundedUpMultipleOf(256) shouldBe MemorySize(2048, MemoryUnit.MB)
    assertThrows[IllegalArgumentException](MemorySize(2000, MemoryUnit.MB).asRoundedUpMultipleOf(0))
    assertThrows[IllegalArgumentException](MemorySize(2000, MemoryUnit.MB).asRoundedUpMultipleOf(-5))
  }
}
