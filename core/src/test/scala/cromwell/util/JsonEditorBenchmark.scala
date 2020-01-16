package cromwell.util

import java.nio.file.{Files, Paths}

import cats.data.NonEmptyList
import cats.syntax.either._
import cromwell.util.ErrorOrUtil._
import cromwell.util.JsonEditor._
import io.circe.Json
import io.circe.parser.parse
import org.scalameter.api._
import org.scalameter.picklers.Implicits._

// Run with:
//   sbt "core/benchmark:testOnly cromwell.util.JsonEditorBenchmark"
// Needs:
//   * A file called bad.json and a file called bec.json in your pwd.
//   * Each json must have at least an "id" field.
object JsonEditorBenchmark extends Bench[Double] {

  /* Mock Json */
  val jsonPool = Map(14 -> "bad.json", 32 -> "bec.json")
  val jsonTestParse = Gen.enumeration("Parse: json size MB")(14, 32)
  val jsonTestExcludeKeys = Gen.enumeration("Exclude keys traversal: json size MB")(14, 32)
  val jsonTestIncludeKeys = Gen.enumeration("Include keys traversal: json size MB")(14, 32)
  val jsonTestAugmentLabels = Gen.enumeration("Augment labels traversal: json size MB")(14, 32)
  val excludeKeys = Map(14 -> "mt_", 32 -> "status")

  def loadJson(fileName: String) : String =  new String(Files.readAllBytes(Paths.get(
    new java.io.File("./engine/src/test/resources")
    .getCanonicalPath, fileName)))

  val jsonStrs = jsonPool map { case (sz, fn) => sz → loadJson(fn) }
  val jsonTree: Map[Int, Either[String, Json]] = jsonStrs map { case (sz, fn) => sz → parse(fn).leftMap(_.toString) }

  /* Benchmark configuration */
  lazy val measurer = new Measurer.Default
  lazy val executor = LocalExecutor(new Executor.Warmer.Default, Aggregator.average, measurer)
  lazy val reporter = new LoggingReporter[Double]
  lazy val persistor = Persistor.None

  // Increase or decrease exec.benchRuns to repeat the same test numerous times.
  // The average measured clock time will be the output.
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
        jsonTree(sz).map { json =>
          updateLabels(json, Map(json.workflowId.get -> Map(("new", "label"))))
        }.right.get // traversal
      }
    }
  }
}

