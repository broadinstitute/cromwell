package cromwell.engine.workflow.tokens.large

import akka.actor.ActorSystem
import akka.testkit.TestProbe
import cromwell.engine.workflow.tokens.{RoundRobinQueueIterator, TokenQueue}
import cromwell.engine.workflow.tokens.TokenQueue.TokenQueuePlaceholder

object TokenDispenserBenchmark {
  implicit val actorSystem: ActorSystem = ActorSystem("TokenDispenserBenchmark")
  val actorToQueue = TestProbe().ref

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

  def useEntireAvailability(tokenQueueIn: TokenQueue, jobsAtATime: Int): TokenQueue = {
    var iterator = new RoundRobinQueueIterator(List(tokenQueueIn), 0, List.empty)

    while (iterator.hasNext) {
      iterator.take(jobsAtATime).toList
      val newQueue = iterator.updatedQueues.head
      iterator = new RoundRobinQueueIterator(List(newQueue), 0, List.empty)
    }

    iterator.updatedQueues.head
  }
}
