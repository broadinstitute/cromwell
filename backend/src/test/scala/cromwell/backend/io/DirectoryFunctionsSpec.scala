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
    override def copyFile(source: String, destination: String) = ???
    override def glob(pattern: String) = ???
    override def size(path: String) = ???
    override def readFile(path: String, maxBytes: Option[Int], failOnOverflow: Boolean)  = ???
    override def pathFunctions = ???
    override def writeFile(path: String, content: String) = ???
    override implicit def ec = ???
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
