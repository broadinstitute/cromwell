package cromwell.binding.values

import cromwell.engine.Hashing._
import cromwell.engine.backend.local.SharedFileSystem
import cromwell.util.HashUtil
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}

class WdlPrimitiveSpec extends FlatSpec with Matchers with HashUtil {

  it should "produce correct Hash" in {
    implicit val fileHasher = SharedFileSystem.sharedFsFileHasher

    // Object with same hashes
    val sameHashTable = Table(
      ("p1", "p2"),

      (WdlString("a"), WdlString("a")),
      (WdlFloat(1.0F), WdlFloat(1.0F)),
      (WdlBoolean(true), WdlBoolean(true)),
      (WdlBoolean(false), WdlBoolean(false)),
      (WdlInteger(1), WdlInteger(1))
    )

    val differentHashTable = Table(
      ("p1", "p2"),

      (WdlString("a"), WdlString("b")),
      (WdlString("1"), WdlFloat(1.0F)),
      (WdlString("1"), WdlInteger(1)),
      (WdlString("true"), WdlBoolean(true)),
      (WdlString("false"), WdlBoolean(false)),
      (WdlBoolean(true), WdlBoolean(false)),
      (WdlInteger(1), WdlInteger(2)),
      (WdlFloat(1.0F), WdlFloat(2.0F)),
      (WdlInteger(1), WdlFloat(1.0F))
    )

    forAll(sameHashTable) { (p1, p2) => p1.computeHash should be(p2.computeHash) }
    forAll(differentHashTable) { (p1, p2) => p1.computeHash shouldNot be(p2.computeHash) }
  }
}
