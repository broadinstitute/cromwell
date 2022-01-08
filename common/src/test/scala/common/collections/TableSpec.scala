package common.collections

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class TableSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "Table"

  val someTable = Table(
    Map(
      "a" -> Map("b" -> "c"),
      "d" -> Map("b" -> "f")
    )
  )

  it should "create an empty Table" in {
    Table.empty.table shouldBe empty
  }

  it should "fill a table" in {
    Table.fill(List(
      ("a", "b", "c"),
      ("d", "e", "f")
    )).table shouldBe Map(
      "a" -> Map("b" -> "c"),
      "d" -> Map("e" -> "f")
    )
  }

  it should "implement contains" in {
    someTable.contains("a", "b") shouldBe true
    someTable.contains("a", "a") shouldBe false
    someTable.contains("g", "a") shouldBe false
    someTable.contains("g", "b") shouldBe false
  }

  it should "implement getValue" in {
    someTable.getValue("a", "b") shouldBe Some("c")
    someTable.getValue("a", "a") shouldBe None
    someTable.getValue("g", "a") shouldBe None
    someTable.getValue("g", "b") shouldBe None
  }

  it should "implement rowOptional" in {
    someTable.rowOptional("a") shouldBe Some(Map("b" -> "c"))
    someTable.rowOptional("b") shouldBe None
  }
  
  it should "implement row" in {
    someTable.row("a") shouldBe Map("b" -> "c")
    someTable.row("b") shouldBe empty
  }

  it should "implement column" in {
    someTable.column("b") shouldBe Map("a" -> "c", "d" -> "f")
    someTable.column("a") shouldBe empty
  }

  it should "implement add" in {
    someTable.add("0", "1", "2") shouldBe Table(
      Map(
        "a" -> Map("b" -> "c"),
        "d" -> Map("b" -> "f"),
        "0" -> Map("1" -> "2")
      )
    )
  }

  it should "implement addAll" in {
    someTable.addAll(
      List(
        ("0", "1", "2"),
        ("3", "2", "4")
      )
    ) shouldBe Table(
      Map(
        "a" -> Map("b" -> "c"),
        "d" -> Map("b" -> "f"),
        "0" -> Map("1" -> "2"),
        "3" -> Map("2" -> "4")
      )
    )
  }

  it should "implement addTriplet" in {
    someTable.addTriplet(
        ("0", "1", "2")
    ) shouldBe Table(
      Map(
        "a" -> Map("b" -> "c"),
        "d" -> Map("b" -> "f"),
        "0" -> Map("1" -> "2")
      )
    )
  }

  it should "implement valuesTriplet" in {
    someTable.valuesTriplet.toList shouldBe List(
      ("a", "b", "c"),
      ("d", "b", "f")
    )
  }
}
