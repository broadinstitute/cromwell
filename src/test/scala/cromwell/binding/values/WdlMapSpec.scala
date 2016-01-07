package cromwell.binding.values

import cromwell.binding.types._
import cromwell.engine.Hashing._
import cromwell.engine.backend.local.SharedFileSystem
import cromwell.util.HashUtil
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}

class WdlMapSpec extends FlatSpec with Matchers with HashUtil {

  it should "produce correct Hash" in {
    implicit val fileHasher = SharedFileSystem.sharedFsFileHasher

    // Note that we don't need to test with every single possible WdlType as WdlMapType, because each type is responsible for producing a valid hash.

    val refMap = WdlMap(WdlMapType(WdlStringType, WdlFileType), Map(string1 -> file1, anotherString -> anotherFile))
    val emptyWdlMap = WdlMap(WdlMapType(WdlStringType, WdlFileType), Map.empty[WdlValue, WdlValue])

    // Object with same hashes
    val sameHashTable = Table(
      ("map1", "map2"),

      // Just in case...
      (refMap, refMap),
      // Empty objects
      (emptyWdlMap, emptyWdlMap),
      // Same inner values but different WdlValues
      (refMap, WdlMap(WdlMapType(WdlStringType, WdlFileType), Map(sameAsString1 -> sameAsfile1, anotherString -> anotherFile))),
      // Shuffled ordering
      (refMap, WdlMap(WdlMapType(WdlStringType, WdlFileType), Map(anotherString -> anotherFile, string1 -> file1)))
    )

    val differentHashTable = Table(
      ("map1", "map2"),

      (refMap, emptyWdlMap),
      // Missing an entry
      (refMap, WdlMap(WdlMapType(WdlStringType, WdlFileType),Map(string1 -> file1))),
      // One more entry
      (refMap, WdlMap(WdlMapType(WdlStringType, WdlFileType),Map(string1 -> file1, anotherString -> anotherFile, WdlString("uninvitedGuest") -> anotherFile))),
      // With a different key
      (refMap, WdlMap(WdlMapType(WdlStringType, WdlFileType), Map(WdlString("whoAmI") -> file1, anotherString -> anotherFile))),
      // With a different value
      (refMap, WdlMap(WdlMapType(WdlStringType, WdlFileType), Map(string1 -> anotherFile, anotherString -> anotherFile)))
    )

    forAll(sameHashTable) { (obj1, obj2) => obj1.computeHash should be(obj2.computeHash) }
    forAll(differentHashTable) { (obj1, obj2) => obj1.computeHash shouldNot be(obj2.computeHash) }
  }

}
