package cromwell.backend.io

import org.scalatest.{FlatSpec, Matchers}
import better.files._
import cromwell.core.path.DefaultPathBuilder
import wom.expression.IoFunctionSet.{IoDirectory, IoFile}

import scala.concurrent.Await
import scala.concurrent.duration.Duration
class DirectoryFunctionsSpec extends FlatSpec with Matchers {
  behavior of "DirectoryFunctions"

  val functions = new DirectoryFunctions {
    override def pathBuilders = List(DefaultPathBuilder)
    override def copyFile(source: String, destination: String) = throw new UnsupportedOperationException()
    override def glob(pattern: String) = throw new UnsupportedOperationException()
    override def size(path: String) = throw new UnsupportedOperationException()
    override def readFile(path: String, maxBytes: Option[Int], failOnOverflow: Boolean)  = throw new UnsupportedOperationException()
    override def pathFunctions = throw new UnsupportedOperationException()
    override def writeFile(path: String, content: String) = throw new UnsupportedOperationException()
    override implicit def ec = throw new UnsupportedOperationException()
    override def createTemporaryDirectory(name: Option[String]) = throw new UnsupportedOperationException()
    override def asyncIo = throw new UnsupportedOperationException()
  }

  "listDirectory" should "exclude visited directories when listing" in {
    val testDir = File.newTemporaryDirectory()
    val rootDir = (testDir / "rootDir").createDirectories()
    val innerDir = (rootDir / "innerDir").createDirectories()
    val link = innerDir / "linkToRootDirInInnerDir"
    link.symbolicLinkTo(rootDir)
    
    def listRecursively(path: String)(visited: Vector[String] = Vector.empty): Iterator[String] = {
      Await.result(functions.listDirectory(path)(visited), Duration.Inf) flatMap {
        case IoFile(v) => List(v)
        case IoDirectory(v) => List(v) ++ listRecursively(v)(visited :+ path)
      }
    }

    listRecursively(rootDir.pathAsString)().toList shouldBe List(innerDir, link).map(_.pathAsString)
  }
}
