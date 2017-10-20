package centaur

import centaur.test.standard.CentaurTestFormat.RestartFormat
import org.scalatest._

import scala.language.postfixOps

@DoNotDiscover
class StandardTestCaseSpec(cromwellBackends: List[String]) extends AbstractCentaurTestCaseSpec(cromwellBackends) with ParallelTestExecution {
  
  def this() = this(CentaurTestSuite.cromwellBackends)
  
  allTestCases.filter( _.testFormat match {
      case _: RestartFormat => false
      case _ => true 
    }
  ) foreach executeStandardTest
}
