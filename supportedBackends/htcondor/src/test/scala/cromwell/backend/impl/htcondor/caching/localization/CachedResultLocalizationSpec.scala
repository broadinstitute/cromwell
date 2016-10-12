package cromwell.backend.impl.htcondor.caching.localization

import java.nio.file.Files

import cromwell.core.{JobOutput, CallOutputs}
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import wdl4s.types.{WdlArrayType, WdlFileType}
import wdl4s.values.{WdlArray, WdlSingleFile, WdlString}

class CachedResultLocalizationSpec extends WordSpecLike with Matchers with BeforeAndAfterAll {
  private class CachedResultLocalizationMock extends CachedResultLocalization
  private val defaultTmpDir = Files.createTempDirectory("cachedFiles").toAbsolutePath
  private val defaultCachedFile = defaultTmpDir.resolve("input.txt")
  private val newTmpDir = Files.createTempDirectory("newFiles").toAbsolutePath
  private val newTmpFile = newTmpDir.resolve(defaultCachedFile.getFileName())
  private val cachedResults = new CachedResultLocalizationMock()
  private val defaultFileArray = Seq("arrInput1.txt", "arrInput2.txt", "arrInput3.txt").map(defaultTmpDir.resolve(_).toAbsolutePath)
  private val newFileArray = Seq("arrInput1.txt", "arrInput2.txt", "arrInput3.txt").map(newTmpDir.resolve(_).toAbsolutePath)

  override def afterAll() = {
    Seq(defaultCachedFile, newTmpFile) ++ newFileArray ++ Seq(defaultTmpDir, newTmpDir) foreach { _.toFile.delete() }
  }

  "CachedResultLocalization" should {
    "localize file path via symbolic link" in {
      val slPath = cachedResults.localizePathViaSymbolicLink(defaultCachedFile, newTmpFile)
      assert(Files.isSymbolicLink(slPath))
      Files.delete(newTmpFile)
    }

    "not localize dir path via symbolic link" in {
      assertThrows[UnsupportedOperationException](cachedResults.localizePathViaSymbolicLink(defaultTmpDir, newTmpFile))
    }

    "localize cached job outputs which are WDL files using symbolic link" in {
      val outputs: CallOutputs = Map("File1" -> JobOutput(WdlSingleFile(defaultCachedFile.toAbsolutePath.toString)))
      val newJobOutputs = cachedResults.localizeCachedOutputs(newTmpDir, outputs)
      newJobOutputs foreach { case (lqn, jobOutput) =>
        assert(jobOutput.wdlValue.valueString == newTmpFile.toString)
      }
    }

    "localize cached job outputs which are WDL File Array using symbolic link" in {
      val wdlArray = WdlArray(WdlArrayType(WdlFileType), defaultFileArray.map(file => WdlSingleFile(file.toString())))
      val outputs = Map("File1" -> JobOutput(wdlArray))
      val newJobOutputs = cachedResults.localizeCachedOutputs(newTmpDir, outputs)
      newJobOutputs foreach { case (lqn, jobOutput) =>
        val wdlArray = jobOutput.wdlValue.asInstanceOf[WdlArray].value
        wdlArray foreach { entry =>
          assert(!entry.valueString.contains(defaultTmpDir.toString))
          assert(entry.valueString.contains(newTmpDir.toString))
        }
      }
    }

    "not localize cached job outputs which are not WDL files" in {
      val outputs = Map("String1" -> JobOutput(WdlString(defaultCachedFile.toAbsolutePath.toString)))
      val newJobOutputs = cachedResults.localizeCachedOutputs(newTmpDir, outputs)
      newJobOutputs foreach { case (lqn, jobOutput) =>
        assert(jobOutput.wdlValue.valueString == defaultCachedFile.toString)
      }
    }
  }
}
