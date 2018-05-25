package wdl.draft3.transforms.wom2wdlom

import cats.implicits._
import java.nio.file.Paths
import org.scalatest.{FlatSpec, Matchers}

import wdl.draft2.model.WdlNamespace
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomBundleMakers.wdlDraft2NamespaceWomBundleMaker

import wdl.draft3.transforms.parsing.{FileStringParserInput, stringToAst}
import wdl.draft3.transforms.wom2wdlom.WomToWdlomImpl.womBundleToFileElement
import wdl.draft3.transforms.wom2wdlom.WomToWdlom.ops._
import wdl.draft3.transforms.wdlom2wdl.WdlWriter.ops._
import wdl.draft3.transforms.wdlom2wdl.WdlWriterImpl.fileElementWriter
import wdl.draft3.transforms.ast2wdlom.astToFileElement
import wdl.draft3.transforms.wdlom2wom.FileElementToWomBundle.fileElementToWomBundle
import wdl.draft3.transforms.wdlom2wom.FileElementToWomBundleInputs

import wom.executable.WomBundle

class WomToWdlomSpec extends FlatSpec with Matchers {

  val testFiles = Seq(
    "scripts/test_upgrade/scatter_files.wdl",
    "centaur/src/main/resources/standardTestCases/taskless_engine_functions/taskless_engine_functions.wdl"
  ).map(Paths.get(_))

  testFiles foreach { file =>
    it should s"upgrade file ${file.getFileName} into something that parses back into WOM" in {
      val wdl: WdlNamespace = WdlNamespace.loadUsingPath(file, None, None).get

      val womBundle: WomBundle = wdlDraft2NamespaceWomBundleMaker.toWomBundle(wdl).right.get

      val sourceString = womBundle.toWdlom.toWdlV1

      val v1objectModel = (stringToAst andThen astToFileElement).run(FileStringParserInput(sourceString, file.getFileName.toString)).right.get

      val generatedWomBundle = fileElementToWomBundle.toWomBundle(FileElementToWomBundleInputs(v1objectModel, "", List.empty, List.empty)).right.get

      generatedWomBundle shouldEqual womBundle
    }
  }

}
