package cromwell.util

import java.nio.file.{Files, Paths}

import JsonEditor._
import cats.data.NonEmptyList
import cats.syntax.either._
import io.circe.parser.parse
import org.scalameter.api._
import org.scalameter.picklers.Implicits._

object JsonEditorBenchmark extends Bench[Double] {

  /* Mock Json */
  val jsonPool = Map(14 -> "bad.json", 32 -> "bec.json")
  val jsonTestParse = Gen.enumeration("Parse: json size MB")(14, 32)
  val jsonTestExcludeKeys = Gen.enumeration("Exclude keys traversal: json size MB")(14, 32)
  val jsonTestIncludeKeys = Gen.enumeration("Include keys traversal: json size MB")(14, 32)
  val jsonTestAugmentLabels = Gen.enumeration("Augment labels traversal: json size MB")(14, 32)
  val excludeKeys = Map(14 -> "mt_", 32 -> "status")

  def loadJson(fileName: String) : String =  new String(Files.readAllBytes(Paths.get(
    new java.io.File(".")
    .getCanonicalPath, fileName)))

  val jsonStrs = jsonPool map { case (sz, fn) => sz → loadJson(fn) }
  val jsonTree = jsonStrs map { case (sz, fn) => sz → parse(fn).leftMap(_.toString) }

  /* Benchmark configuration */
  lazy val measurer = new Measurer.Default
  lazy val executor = LocalExecutor(new Executor.Warmer.Default, Aggregator.average, measurer)
  lazy val reporter = new LoggingReporter[Double]
  lazy val persistor = Persistor.None

  // Increase or decrease exec.benchRuns to repeat the same test numerous times.
  // The avaerage measured clock time will be the output.
  performance of "JsonEditor with circe" in {
    // Measures how fast the JsonEditor can exclude keys using circe.
    measure method "parse" in {
      using(jsonTestParse) config (
        exec.maxWarmupRuns -> 1,
        exec.benchRuns -> 1
      ) in { sz =>
        for {
          tree <- parse(jsonStrs(sz)).leftMap(_.toString)
        } yield tree
      }
    }

    measure method "includeExcludeJson(_, None, Some(NonEmptyList.of(<some exclude key>)))" in {
      using(jsonTestExcludeKeys) config (
        exec.maxWarmupRuns -> 1,
        exec.benchRuns -> 1
      ) in { sz =>
        for {
          tree <- jsonTree(sz) // already parsed
          newTree = includeExcludeJson(tree, None, Some(NonEmptyList.of(excludeKeys(sz)))) // traversal
        } yield newTree
      }
    }

    measure method "includeJson(_, NonEmptyList.one(\"message\"))" in {
      using(jsonTestIncludeKeys) config(
        exec.maxWarmupRuns -> 1,
        exec.benchRuns -> 1
      ) in { sz =>
        jsonTree(sz).map(includeJson(_, NonEmptyList.one("message"))).right.get // traversal
      }
    }

    measure method "augmentLabels(_, Map((\"new\",\"label\")))" in {
      using(jsonTestAugmentLabels) config (
        exec.maxWarmupRuns -> 1,
        exec.benchRuns -> 1
      ) in { sz =>
        jsonTree(sz).map(augmentLabels(_, Map(("new","label")))).right.get // traversal
      }
    }
  }
}

