package centaur.reporting

/**
  * Information about a test.
  *
  * @param testName The test name.
  * @param retries The total number of retries.
  * @param attempt The zero based attempt.
  */
case class TestEnvironment(testName: String, retries: Int, attempt: Int)
