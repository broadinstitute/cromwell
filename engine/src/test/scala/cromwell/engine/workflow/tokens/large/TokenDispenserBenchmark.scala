package cromwell.engine.workflow.tokens.large

import akka.actor.ActorSystem
import akka.testkit.TestProbe
import cromwell.core.ExecutionStatus.{apply => _}
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.{NullTokenEventLogger, RoundRobinQueueIterator, TokenQueue}
import cromwell.engine.workflow.tokens.TokenQueue.TokenQueuePlaceholder
import org.scalameter.api._
import org.scalameter.picklers.Implicits._
import spray.json.DefaultJsonProtocol

/**
  * Benchmarks the performance of the execution store using ScalaMeter (http://scalameter.github.io/)
  * This is not run automatically by "sbt test". To run this test specifically, either use intellij integration, or run
  * sbt "engine/benchmark:testOnly cromwell.engine.workflow.tokens.large.TokenDispenserBenchmark"
  * sbt benchmark:test will run all ScalaMeter tests
  */
object TokenDispenserBenchmark extends Bench[Double] with DefaultJsonProtocol {
  implicit val actorSystem: ActorSystem = ActorSystem("TokenDispenserBenchmark")
  val actorToQueue = TestProbe().ref
  val ScaleFactor = 10000

  /* Benchmark configuration */
  override lazy val measurer = new Measurer.Default
//  override lazy val executor = SeparateJvmsExecutor(new Executor.Warmer.Default, Aggregator.average, measurer)
  override lazy val executor = LocalExecutor(new Executor.Warmer.Default, Aggregator.average, measurer)
  override lazy val reporter = new LoggingReporter[Double]
  override lazy val persistor = Persistor.None

  def fillQueue(tokenQueueIn: TokenQueue, jobsPerGroup: Int, hogGroups: Int): TokenQueue = {
    val hogGroupNames = (0 until hogGroups) map { i => s"hogGroup$i" }
    fillQueue(tokenQueueIn, jobsPerGroup, hogGroupNames.toList)
  }

  def fillQueue(tokenQueueIn: TokenQueue, jobsPerGroup: Int, hogGroups: List[String]): TokenQueue = {
    var tokenQueue = tokenQueueIn
    hogGroups foreach { hogGroup =>
    (0 until jobsPerGroup) foreach { _ =>
        tokenQueue = tokenQueue.enqueue(TokenQueuePlaceholder(actorToQueue, hogGroup))
      }
    }
    tokenQueue
  }

  def useQueue(tokenQueueIn: TokenQueue, jobsToDequeue: Int, jobsAtATime: Int) = {
    var tokenQueue = tokenQueueIn
    (0 until (jobsToDequeue / jobsAtATime)) foreach { _ =>
      val iterator = new RoundRobinQueueIterator(List(tokenQueue), 0)
      iterator.take(jobsAtATime)
      tokenQueue = iterator.updatedQueues.head
    }
    tokenQueue
  }

  def useEntireAvailability(tokenQueueIn: TokenQueue, jobsAtATime: Int): TokenQueue = {
    var iterator = new RoundRobinQueueIterator(List(tokenQueueIn), 0)

    while (iterator.hasNext) {
      iterator.take(jobsAtATime).toList
      val newQueue = iterator.updatedQueues.head
      iterator = new RoundRobinQueueIterator(List(newQueue), 0)
    }

    iterator.updatedQueues.head
  }

  performance of "RoundRobinQueueIterator" in {

    measure method "enqueuing and dequeuing with multiple hog groups" in {
      val poolSize = 5 * ScaleFactor
      val jobCounts: Gen[Int] = Gen.range("initialJobsInQueue")(from = 1 * ScaleFactor, upto = 15 * ScaleFactor, hop = 3 * ScaleFactor)
      val jobsAtATime = 50

      val queues = for {
        jobCount <- jobCounts
        hogFactor <- Gen.exponential("hogFactor")(from = 1, until = Math.min(jobCount, 65536), factor = 4)
        tokenType = JobExecutionTokenType("BIG_PAPI", Some(poolSize), hogFactor)

      } yield (jobCount, hogFactor, tokenType)

      using(queues) in { case (jobCount, hogFactor, tokenType) =>
        var tokenQueue = TokenQueue(tokenType, NullTokenEventLogger)
        val hogGroups = ((0 until hogFactor) map { i => i -> s"hogGroup$i" }).toMap

        tokenQueue = fillQueue(tokenQueue, jobCount / hogFactor, hogGroups.values.toList)
        useEntireAvailability(tokenQueue, jobsAtATime)
      }

    }

  }
}
