package centaur.reporting

import centaur.test.standard.CentaurTestCase

/**
  * Information about a test.
  *
  * @param testCase The Centaur test case.
  * @param retries The total number of retries.
  * @param attempt The zero based attempt.
  */
case class TestEnvironment(testCase: CentaurTestCase, retries: Int, attempt: Int)
