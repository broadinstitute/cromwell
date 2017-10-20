package cromwell.backend

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers, TryValues}
import wdl4s.parser.MemoryUnit

class MemorySizeSpec extends FlatSpec with Matchers with TryValues {
  "MemorySize" should "stringify properly for integer values" in {
    val memTable = Table(
      ("memorySize", "memorySizeString"),
      (MemorySize(10, MemoryUnit.Bytes), "10 B"),
      (MemorySize(10, MemoryUnit.KB), "10 KB"),
      (MemorySize(10, MemoryUnit.MB), "10 MB"),
      (MemorySize(10, MemoryUnit.GB), "10 GB"),
      (MemorySize(10, MemoryUnit.TB), "10 TB"),
      (MemorySize(10, MemoryUnit.KiB), "10 KiB"),
      (MemorySize(10, MemoryUnit.MiB), "10 MiB"),
      (MemorySize(10, MemoryUnit.GiB), "10 GiB"),
      (MemorySize(10, MemoryUnit.TiB), "10 TiB")
    )

    forAll(memTable) { (memorySize, memorySizeString) =>
      memorySize.toString shouldEqual memorySizeString
    }
  }

  it should "stringify properly for non-integer values" in {
    val memTable = Table(
      ("memorySize", "memorySizeString"),
      (MemorySize(10.5, MemoryUnit.KB), "10.5 KB"),
      (MemorySize(10.5, MemoryUnit.MB), "10.5 MB"),
      (MemorySize(10.5, MemoryUnit.GB), "10.5 GB"),
      (MemorySize(10.5, MemoryUnit.TB), "10.5 TB"),
      (MemorySize(10.5, MemoryUnit.KiB), "10.5 KiB"),
      (MemorySize(10.5, MemoryUnit.MiB), "10.5 MiB"),
      (MemorySize(10.5, MemoryUnit.GiB), "10.5 GiB"),
      (MemorySize(10.5, MemoryUnit.TiB), "10.5 TiB")
    )

    forAll(memTable) { (memorySize, memorySizeString) =>
      memorySize.toString shouldEqual memorySizeString
    }
  }

  it should "convert to bytes properly" in {
    val memTable = Table(
      ("memorySize", "bytes"),
      (MemorySize(10.5, MemoryUnit.KB), 10500.0),
      (MemorySize(10.5, MemoryUnit.MB), 10500000.0),
      (MemorySize(10.5, MemoryUnit.GB), 10500000000.0),
      (MemorySize(10.5, MemoryUnit.TB), 10500000000000.0),
      (MemorySize(10.5, MemoryUnit.KiB), 10752.0),
      (MemorySize(10.5, MemoryUnit.MiB), 11010048.0),
      (MemorySize(10.5, MemoryUnit.GiB), 11274289152.0),
      (MemorySize(10.5, MemoryUnit.TiB), 11544872091648.0)
    )

    forAll(memTable) { (memorySize, bytes) =>
      memorySize.bytes shouldEqual bytes
    }
  }

  it should "convert to other units properly" in {
    val memTable = Table(
      ("memorySize", "newUnit", "result"),
      (MemorySize(1000, MemoryUnit.Bytes), MemoryUnit.KB, MemorySize(1, MemoryUnit.KB)),
      (MemorySize(1000000, MemoryUnit.Bytes), MemoryUnit.MB, MemorySize(1, MemoryUnit.MB)),
      (MemorySize(1024, MemoryUnit.Bytes), MemoryUnit.KiB, MemorySize(1, MemoryUnit.KiB)),
      (MemorySize(1, MemoryUnit.KB), MemoryUnit.Bytes, MemorySize(1000, MemoryUnit.Bytes)),
      (MemorySize(1, MemoryUnit.MB), MemoryUnit.Bytes, MemorySize(1000000, MemoryUnit.Bytes)),
      (MemorySize(1, MemoryUnit.KiB), MemoryUnit.Bytes, MemorySize(1024, MemoryUnit.Bytes))
    )

    forAll(memTable) { (memorySize, newUnit, result) =>
      memorySize.to(newUnit) shouldEqual result
    }
  }

  it should "parse strings" in {
    val memTable = Table(
      ("unparsed", "memorySize"),
      ("1000 B", MemorySize(1000, MemoryUnit.Bytes)),
      ("100KB", MemorySize(100, MemoryUnit.KB)),
      ("10.2MB", MemorySize(10.2, MemoryUnit.MB)),
      ("10.2MiB", MemorySize(10.2, MemoryUnit.MiB)),
      ("1.5 GB", MemorySize(1.5, MemoryUnit.GB))
    )

    forAll(memTable) { (unparsed, memorySize) =>
      MemorySize.parse(unparsed).get shouldEqual memorySize
    }
  }

  it should "stringify values" in {
    val memTable = Table(
      ("memorySize", "string"),
      (MemorySize(1000, MemoryUnit.Bytes), "1000 B"),
      (MemorySize(100, MemoryUnit.KB), "100 KB"),
      (MemorySize(10.2, MemoryUnit.MB), "10.2 MB"),
      (MemorySize(10.2, MemoryUnit.MiB), "10.2 MiB"),
      (MemorySize(1.5, MemoryUnit.GB), "1.5 GB")
    )

    forAll(memTable) { (memorySize, string) =>
      memorySize.toString shouldEqual string
      MemorySize.parse(string).get shouldEqual memorySize
    }
  }

  it should "parse a unit from its suffix" in {
    val memTable = Table(
      ("memorySize", "string"),
      (MemoryUnit.Bytes, "B"),
      (MemoryUnit.KB, "K"),
      (MemoryUnit.KB, "KB"),
      (MemoryUnit.MB, "M"),
      (MemoryUnit.MB, "MB"),
      (MemoryUnit.GB, "G"),
      (MemoryUnit.GB, "GB"),
      (MemoryUnit.TB, "T"),
      (MemoryUnit.TB, "TB"),
      (MemoryUnit.KiB, "Ki"),
      (MemoryUnit.KiB, "KiB"),
      (MemoryUnit.MiB, "Mi"),
      (MemoryUnit.MiB, "MiB"),
      (MemoryUnit.GiB, "Gi"),
      (MemoryUnit.GiB, "GiB"),
      (MemoryUnit.TiB, "Ti"),
      (MemoryUnit.TiB, "TiB")
    )

    forAll(memTable) { (unit, suffix) =>
      MemoryUnit.fromSuffix(suffix) shouldEqual unit
    }

    an[IllegalArgumentException] should be thrownBy {
      MemoryUnit.fromSuffix("Unknown suffix")
    }
  }
}
