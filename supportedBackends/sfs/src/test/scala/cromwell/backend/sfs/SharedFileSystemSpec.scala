package cromwell.backend.sfs

import akka.actor.ActorContext
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.BackendSpec
import cromwell.core.CromwellFatalExceptionMarker
import cromwell.core.path.{DefaultPathBuilder, Path}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wom.values.WomSingleFile

import scala.io.Source

class SharedFileSystemSpec extends FlatSpec with Matchers with Mockito with TableDrivenPropertyChecks with BackendSpec {

  behavior of "SharedFileSystem"

  val defaultLocalization = ConfigFactory.parseString(""" localization: [copy, hard-link, soft-link] """)
  val hardLinkLocalization = ConfigFactory.parseString(""" localization: [hard-link] """)
  val softLinkLocalization = ConfigFactory.parseString(""" localization: [soft-link] """)
  val cachedCopyLocalization = ConfigFactory.parseString(""" localization: [cached-copy] """)
  val cachedCopyLocalizationMaxHardlinks = ConfigFactory.parseString("""{localization: [cached-copy], max-hardlinks: 3 }""")
  val localPathBuilder = List(DefaultPathBuilder)


  def localizationTest(config: Config,
                       docker: Boolean,
                       fileInCallDir: Boolean = false,
                       fileAlreadyExists: Boolean = false,
                       symlink: Boolean = false,
                       cachedCopy: Boolean = false,
                       linkNb: Int = 1) = {
    val callDir = DefaultPathBuilder.createTempDirectory("SharedFileSystem")
    val orig = if (fileInCallDir) callDir.createChild("inputFile") else DefaultPathBuilder.createTempFile("inputFile")
    val dest = if (fileInCallDir) orig else callDir./(orig.parent.pathAsString.hashCode.toString())./(orig.name)
    val testText =
      """This is a simple text to check if the localization
        | works correctly for the file contents.
        |""".stripMargin
    orig.touch()
    orig.writeText(testText)
    if (fileAlreadyExists) {
      dest.parent.createPermissionedDirectories()
      dest.touch()
      dest.writeText(testText)
    }

    val inputs = fqnWdlMapToDeclarationMap(Map("input" -> WomSingleFile(orig.pathAsString)))
    val sharedFS = new SharedFileSystem {
      override val pathBuilders = localPathBuilder
      override val sharedFileSystemConfig = config
      override implicit def actorContext: ActorContext = null
      override lazy val cachedCopyDir = Some(DefaultPathBuilder.createTempDirectory("cached-copy"))
    }
    val cachedFile: Option[Path] = sharedFS.cachedCopyDir.map(
      _./(orig.parent.pathAsString.hashCode.toString)./(orig.lastModifiedTime.toEpochMilli.toString + orig.name))
    val localizedinputs = Map(inputs.head._1 -> WomSingleFile(dest.pathAsString))
    val result = sharedFS.localizeInputs(callDir, docker = docker)(inputs)

    result.isSuccess shouldBe true
    result.get.toList should contain theSameElementsAs localizedinputs

    dest.exists shouldBe true
    Source.fromFile(dest.toFile).mkString shouldBe testText
    countLinks(dest) should be(linkNb)
    isSymLink(dest) should be(symlink)

    cachedFile.foreach(_.exists should be(cachedCopy))
    orig.delete(swallowIOExceptions = true)
    dest.delete(swallowIOExceptions = true)
  }

  it should "not localize a file already in the call root" in {
    localizationTest(defaultLocalization, docker = false, fileInCallDir = true)
    localizationTest(defaultLocalization, docker = true, fileInCallDir = true)
  }

  it should "not localize a file already localized" in {
    localizationTest(defaultLocalization, docker = false, fileAlreadyExists = true)
    localizationTest(defaultLocalization, docker = true, fileAlreadyExists = true)
  }

  it should "localize a file via copy" in {
    localizationTest(defaultLocalization, docker = false)
    localizationTest(defaultLocalization, docker = true)
  }

  it should "localize a file via hard link" in {
    localizationTest(hardLinkLocalization, docker = false, linkNb = 2)
    localizationTest(hardLinkLocalization, docker = true, linkNb = 2)
  }

  it should "localize a file via symbolic link" in {
    localizationTest(softLinkLocalization, docker = false, symlink = true)
  }

  it should "localize a file via cached copy" in {
    localizationTest(cachedCopyLocalization, docker = false, cachedCopy = true, linkNb = 2)
    localizationTest(cachedCopyLocalization, docker = true, cachedCopy = true, linkNb = 2)
  }

  it should "throw a fatal exception if localization fails" in {
    val callDir = DefaultPathBuilder.createTempDirectory("SharedFileSystem")
    val orig = DefaultPathBuilder.get("/made/up/origin")

    val inputs = fqnWdlMapToDeclarationMap(Map("input" -> WomSingleFile(orig.pathAsString)))
    val sharedFS = new SharedFileSystem {
      override val pathBuilders = localPathBuilder
      override val sharedFileSystemConfig = defaultLocalization
      override implicit def actorContext: ActorContext = null
    }
    val result = sharedFS.localizeInputs(callDir, docker = false)(inputs)
    result.isFailure shouldBe true
    result.failed.get.isInstanceOf[CromwellFatalExceptionMarker] shouldBe true
  }

  it should "cache only one file if copied multiple times via cached copy" in {
    val callDirs: List[Path] = List(1,2,3).map(x => DefaultPathBuilder.createTempDirectory("SharedFileSystem"))
    val orig = DefaultPathBuilder.createTempFile("inputFile")
    val dests = callDirs.map(_./(orig.parent.pathAsString.hashCode.toString())./(orig.name))
    orig.touch()
    val inputs = fqnWdlMapToDeclarationMap(Map("input" -> WomSingleFile(orig.pathAsString)))
    val sharedFS = new SharedFileSystem {
      override val pathBuilders = localPathBuilder
      override val sharedFileSystemConfig = cachedCopyLocalization
      override implicit def actorContext: ActorContext = null
      override lazy val cachedCopyDir = Some(DefaultPathBuilder.createTempDirectory("cached-copy"))
    }
    val cachedFile: Option[Path] = sharedFS.cachedCopyDir.map(
      _./(orig.parent.pathAsString.hashCode.toString)./(orig.lastModifiedTime.toEpochMilli.toString + orig.name))

    val results = callDirs.map(sharedFS.localizeInputs(_, docker = true)(inputs))

    results.foreach(_.isSuccess shouldBe true)
    dests.foreach(_.exists shouldBe true)
    dests.foreach(countLinks(_) shouldBe 4)

    cachedFile.foreach(_.exists shouldBe true)
    cachedFile.foreach(countLinks(_) shouldBe 4)
    orig.delete(swallowIOExceptions = true)
    dests.foreach(_.delete(swallowIOExceptions = true))
  }

it should "copy the file again when the copy-cached file has exceeded the maximum number of hardlinks" in {
    val callDirs: IndexedSeq[Path] = 1 to 3 map { _ => DefaultPathBuilder.createTempDirectory("SharedFileSystem") }
    val orig = DefaultPathBuilder.createTempFile("inputFile")
    val dests = callDirs.map(_./(orig.parent.pathAsString.hashCode.toString())./(orig.name))
    orig.touch()
    val inputs = fqnWdlMapToDeclarationMap(Map("input" -> WomSingleFile(orig.pathAsString)))
    val sharedFS = new SharedFileSystem {
      override val pathBuilders = localPathBuilder
      override val sharedFileSystemConfig = cachedCopyLocalizationMaxHardlinks
      override implicit def actorContext: ActorContext = null
      override lazy val cachedCopyDir = Some(DefaultPathBuilder.createTempDirectory("cached-copy"))
    }
    val cachedFile: Option[Path] = sharedFS.cachedCopyDir.map(
      _./(orig.parent.pathAsString.hashCode.toString)./(orig.lastModifiedTime.toEpochMilli.toString + orig.name))

    val results = callDirs.map(sharedFS.localizeInputs(_, docker = true)(inputs))

    results.foreach(_.isSuccess shouldBe true)
    dests.foreach(_.exists shouldBe true)
    dests.foreach(countLinks(_) should be <= 3)

    cachedFile.foreach(_.exists shouldBe true)
    cachedFile.foreach(countLinks(_) should be <= 3)
    orig.delete(swallowIOExceptions = true)
    dests.foreach(_.delete(swallowIOExceptions = true))
  }


  private[this] def countLinks(file: Path): Int = file.getAttribute("unix:nlink").asInstanceOf[Int]

  private[this] def isSymLink(file: Path): Boolean = file.isSymbolicLink
}
