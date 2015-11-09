package cromwell.binding.values

import cromwell.binding.types.{WdlStringType, WdlMapType, WdlFileType, WdlArrayType}
import cromwell.engine.backend.local.SharedFileSystem
import cromwell.util.HashUtil
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{Matchers, FlatSpec}

class WdlArraySpec extends FlatSpec with Matchers with HashUtil {

  it should "produce correct Hash" in {
    implicit val fileHasher = SharedFileSystem.sharedFSFileHasher

    // Note that we don't need to test with every single possible WdlType as WdlArrayType, because each type is responsible for producing a valid hash.

    val refArray = WdlArray(WdlArrayType(WdlFileType), Seq(file1, anotherFile))
    val emptyArray = WdlArray(WdlArrayType(WdlFileType), Seq.empty)

    // Object with same hashes
    val sameHashTable = Table(
      ("arr1", "arr2"),

      // Just in case...
      (refArray, refArray),
      // Empty arrays
      (emptyArray, emptyArray),
      // Same inner values but different WdlValues
      (refArray, WdlArray(WdlArrayType(WdlFileType), Seq(sameAsfile1, anotherFile)))
    )

    val differentHashTable = Table(
      ("arr1", "arr2"),

      (refArray, emptyArray),
      // Missing an entry
      (refArray, WdlArray(WdlArrayType(WdlFileType), Seq(file1))),
      // One more entry
      (refArray, WdlArray(WdlArrayType(WdlFileType), Seq(file1, file1, anotherFile))),
      // With a different order
      (refArray, WdlArray(WdlArrayType(WdlFileType), Seq(anotherFile, file1))),
      // A WdlMap with same values
      (refArray, WdlMap(WdlMapType(WdlFileType, WdlFileType), Map(file1 -> anotherFile)))
    )

    forAll(sameHashTable) { (arr1, arr2) => arr1.getHash should be(arr2.getHash) }
    forAll(differentHashTable) { (arr1, arr2) => arr1.getHash shouldNot be(arr2.getHash) }
  }

}
