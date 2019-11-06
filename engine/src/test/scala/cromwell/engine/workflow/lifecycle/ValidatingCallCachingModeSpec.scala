package cromwell.engine.workflow.lifecycle

import cats.data.Validated.{Invalid, Valid}
import cromwell.core.WorkflowOptions
import cromwell.core.callcaching.CallCachingMode
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{Assertion, FlatSpec, Matchers}

import scala.util.{Success, Try}

class ValidatingCallCachingModeSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  def makeOptions(writeOpt: Option[Boolean], readOpt: Option[Boolean]) = {
    val writeValue = writeOpt map { v => s""""write_to_cache": $v""" }
    val readValue = readOpt map { v => s""""read_from_cache": $v""" }
    val workflowOptions = Seq(writeValue, readValue).flatten.mkString(",\n")

    WorkflowOptions.fromJsonString(s"{$workflowOptions}")
  }

  val options = Seq(None, Some(true), Some(false))
  val callCachingEnabled = true
  val invalidBadCaсheResults = true

  val allCombinations = (for {
    writeOption <- options
    readOption <- options
  } yield (makeOptions(writeOption, readOption))).toSet

  // writeCache is ON when config is ON and write_to_cache is None or true
  val writeCacheOnCombinations = (for {
    writeOption <- options if writeOption.isEmpty || writeOption.get
    readOption <- options
  } yield (makeOptions(writeOption, readOption))).toSet

  // readCache is ON when config is ON and read_from_cache is None or true
  val readCacheOnCombinations = (for {
    writeOption <- options
    readOption <- options if readOption.isEmpty || readOption.get
  } yield (makeOptions(writeOption, readOption))).toSet

  val writeCacheOffCombinations = allCombinations -- writeCacheOnCombinations
  val readCacheOffCombinations = allCombinations -- readCacheOnCombinations

  validateCallCachingMode(
    "write cache on options",
    writeCacheOnCombinations,
    callCachingEnabled,
    invalidBadCaсheResults) { _.writeToCache should be(true) }
  validateCallCachingMode(
    "read cache on options",
    readCacheOnCombinations,
    callCachingEnabled,
    invalidBadCaсheResults) { _.readFromCache should be(true) }
  validateCallCachingMode(
    "write cache off options",
    writeCacheOffCombinations,
    callCachingEnabled,
    invalidBadCaсheResults) { _.writeToCache should be(false) }
  validateCallCachingMode(
    "read cache off options",
    readCacheOffCombinations,
    callCachingEnabled,
    invalidBadCaсheResults) { _.readFromCache should be(false) }

  private def validateCallCachingMode(testName: String,
                                      wfOptions: Set[Try[WorkflowOptions]],
                                      callCachingEnabled: Boolean,
                                      invalidBadCacheResults: Boolean)
                                      (verificationFunction: CallCachingMode => Assertion): Unit = {
    it should s"correctly identify $testName" in {
      wfOptions foreach  {
        case Success(wfOptions) =>
          MaterializeWorkflowDescriptorActor.validateCallCachingMode(wfOptions,
            callCachingEnabled, invalidBadCacheResults) match {
            case Valid(activity) => verificationFunction(activity)
            case Invalid(errors) =>
              val errorsList = errors.toList.mkString(", ")
              fail(s"Failure generating Call Config Mode: $errorsList")
          }
        case x => fail(s"Unexpected test tuple: $x")
      }
    }
  }
}
