package cromwell.languages.util

import common.assertion.ErrorOrAssertions._
import cromwell.languages.util.ImportResolver.HttpResolver
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

  it should "resolve an abolute path from a different initial root" in {
    val resolver = HttpResolver(relativeTo = Some("http://abc.com:8000/blah1/blah2/"))
    val toResolve = resolver.pathToLookup("http://def.org:8080/blah3.wdl")
    toResolve shouldBeValid "http://def.org:8080/blah3.wdl"
  }

  it should "resolve a relative path" in {
    val resolver = HttpResolver(relativeTo = Some("http://abc.com:8000/blah1/blah2/"))
    val toResolve = resolver.pathToLookup("tools/cool_tool.wdl")
    toResolve shouldBeValid "http://abc.com:8000/blah1/blah2/tools/cool_tool.wdl"
  }

  it should "resolve a backtracking relative path" in {
    val resolver = HttpResolver(relativeTo = Some("http://abc.com:8000/blah1/blah2/"))
    val toResolve = resolver.pathToLookup("../tools/cool_tool.wdl")
    toResolve shouldBeValid "http://abc.com:8000/blah1/tools/cool_tool.wdl"
  }

  it should "resolve from '/'" in {
    val resolver = HttpResolver(relativeTo = Some("http://abc.com:8000/blah1/blah2/"))
    val toResolve = resolver.pathToLookup("/tools/cool_tool.wdl")
    toResolve shouldBeValid "http://abc.com:8000/tools/cool_tool.wdl"
  }
}
