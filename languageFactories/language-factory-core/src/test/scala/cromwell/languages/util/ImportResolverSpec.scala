package cromwell.languages.util

import java.nio.file.Paths

import common.assertion.ErrorOrAssertions._
import cromwell.core.path.DefaultPath
import cromwell.languages.util.ImportResolver.{DirectoryResolver, HttpResolver}
import org.scalatest.{FlatSpec, Matchers}

class ImportResolverSpec extends FlatSpec with Matchers {
  behavior of "HttpResolver"

  val canon = Map(
    "http://abc.com:8000/blah?x=5&y=10" -> "http://abc.com:8000/blah?x=5&y=10",
    "http://abc.com:8000/bob/loblaw/law/blog/../../blah/blah" -> "http://abc.com:8000/bob/loblaw/blah/blah",
    "http://abc.com:8000/bob/loblaw/law/blog/..//../blah/blah" -> "http://abc.com:8000/bob/loblaw/blah/blah"
  )

  canon foreach { case (from, to) =>
    it should s"canonicalize '$from' correctly as $to" in {
      HttpResolver.canonicalize(from) should be(to)
    }
  }

  val roots = Map(
    "http://abc.com:8000/blah?x=5&y=10" -> "http://abc.com:8000",
    "https://raw.githubusercontent.com/broadinstitute/cromwell/1c41a263e9f251c9fedfec6e6e8cbce2dc3e61c6/centaur/src/main/resources/standardTestCases/biscayne_http_relative_imports.test" -> "https://raw.githubusercontent.com"
  )

  roots foreach { case (from, to) =>
    it should s"root '$from' correctly as $to" in {
      HttpResolver.root(from) should be(to)
    }
  }

  it should "resolve a path from no initial root" in {
    val resolver = HttpResolver()
    val toResolve = resolver.pathToLookup("http://abc.com:8000/blah1/blah2.wdl")
    toResolve shouldBeValid "http://abc.com:8000/blah1/blah2.wdl"
  }

  behavior of "HttpResolver with a 'relativeTo' value"

  val relativeHttpResolver = HttpResolver(relativeTo = Some("http://abc.com:8000/blah1/blah2/"))

  it should "resolve an abolute path from a different initial root" in {
    val pathToLookup = relativeHttpResolver.pathToLookup("http://def.org:8080/blah3.wdl")
    pathToLookup shouldBeValid "http://def.org:8080/blah3.wdl"
  }

  it should "resolve a relative path" in {
    val pathToLookup = relativeHttpResolver.pathToLookup("tools/cool_tool.wdl")
    pathToLookup shouldBeValid "http://abc.com:8000/blah1/blah2/tools/cool_tool.wdl"
  }

  it should "resolve a backtracking relative path" in {
    val pathToLookup = relativeHttpResolver.pathToLookup("../tools/cool_tool.wdl")
    pathToLookup shouldBeValid "http://abc.com:8000/blah1/tools/cool_tool.wdl"
  }

  it should "resolve from '/'" in {
    val pathToLookup = relativeHttpResolver.pathToLookup("/tools/cool_tool.wdl")
    pathToLookup shouldBeValid "http://abc.com:8000/tools/cool_tool.wdl"
  }

  behavior of "directory resolver from root"

  val rootDirectoryResolver = DirectoryResolver(DefaultPath(Paths.get("/")), customName = None)

  it should "resolve a random path" in {
    val pathToLookup = rootDirectoryResolver.resolveAndMakeAbsolute("/path/to/file.wdl")
    pathToLookup shouldBeValid Paths.get("/path/to/file.wdl")
  }

  behavior of "unprotected relative directory resolver"

  val relativeDirectoryResolver = DirectoryResolver(DefaultPath(Paths.get("/path/to/imports/")), customName = None)

  it should "resolve an absolute path" in {
    val pathToLookup = relativeDirectoryResolver.resolveAndMakeAbsolute("/path/to/file.wdl")
    pathToLookup shouldBeValid Paths.get("/path/to/file.wdl")
  }

  it should "resolve a relative path" in {
    val pathToLookup = relativeDirectoryResolver.resolveAndMakeAbsolute("path/to/file.wdl")
    pathToLookup shouldBeValid Paths.get("/path/to/imports/path/to/file.wdl")
  }

  behavior of "protected relative directory resolver"

  val protectedRelativeDirectoryResolver = DirectoryResolver(DefaultPath(Paths.get("/path/to/imports/")), Some("/path/to/imports/"), customName = None)

  it should "resolve a good relative path" in {
    val pathToLookup = protectedRelativeDirectoryResolver.resolveAndMakeAbsolute("path/to/file.wdl")
    pathToLookup shouldBeValid Paths.get("/path/to/imports/path/to/file.wdl")
  }

  it should "not resolve an absolute path" in {
    val pathToLookup = protectedRelativeDirectoryResolver.resolveAndMakeAbsolute("/path/to/file.wdl")
    pathToLookup shouldBeInvalid "/path/to/file.wdl is not an allowed import path"
  }

  it should "not resolve outside of its protected directory" in {
    val pathToLookup = protectedRelativeDirectoryResolver.resolveAndMakeAbsolute("../file.wdl")
    pathToLookup shouldBeInvalid "../file.wdl is not an allowed import path"
  }

}
