package cromwell.services.metadata

import java.util.UUID

import common.assertion.ManyTimes.intWithTimes
import cromwell.core.WorkflowId
import cromwell.services.metadata.impl.ActiveMetadataStatisticsRecorder
import cromwell.services.metadata.impl.MetadataStatisticsRecorder.HeavyMetadataAlert
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wom.values.{WomInteger, WomString}

import scala.util.Random

class MetadataStatisticsRecorderSpec extends AnyFlatSpec with Matchers {

  behavior of "MetadataStatisticsRecorder"

  def uninterestingWriteEvent(workflowId: WorkflowId)(): MetadataEvent = MetadataEvent(MetadataKey(workflowId, None, "foo"), MetadataValue(WomInteger(1)))
  def parentNotificationEvent(rootWorkflowId: WorkflowId, parentWorkflowId: WorkflowId, subWorkflowId: WorkflowId): Vector[MetadataEvent] = Vector(
    MetadataEvent(MetadataKey(subWorkflowId, None, "rootWorkflowId"), MetadataValue(WomString(rootWorkflowId.toString))),
    MetadataEvent(MetadataKey(subWorkflowId, None, "parentWorkflowId"), MetadataValue(WomString(parentWorkflowId.toString)))
  )

  it should "count rows for one workflow and create alerts every interval" in {
    val recorder = new ActiveMetadataStatisticsRecorder(10, 10)
    val workflowId = WorkflowId(UUID.randomUUID())

    (1 to 10) foreach { i =>
      recorder.processEventsAndGenerateAlerts(9 of uninterestingWriteEvent(workflowId)) should be(Vector.empty)
      recorder.processEventsAndGenerateAlerts(1 of uninterestingWriteEvent(workflowId)) should be(Vector(HeavyMetadataAlert(workflowId, 10L * i)))
      ()
    }
  }

  it should "be able to reset counters from odd values" in {
    val recorder = new ActiveMetadataStatisticsRecorder(10, 10)
    val workflowId = WorkflowId(UUID.randomUUID())

    var runningCounter: Long = 0L
    (1 to 10) foreach { i =>
      val bigMetadataDump = Math.abs(Random.nextInt(101)) + 10
      val expectedCountAfterProcessing = runningCounter + bigMetadataDump

      recorder.processEventsAndGenerateAlerts(bigMetadataDump of uninterestingWriteEvent(workflowId)) should be(Vector(HeavyMetadataAlert(workflowId, expectedCountAfterProcessing)))
      runningCounter = expectedCountAfterProcessing
      ()
    }
  }

  it should "count rows for multiple workflow and create alerts every interval" in {
    val recorder = new ActiveMetadataStatisticsRecorder(10, 10)
    val workflowId1 = WorkflowId(UUID.randomUUID())
    val workflowId2 = WorkflowId(UUID.randomUUID())
    val workflowId3 = WorkflowId(UUID.randomUUID())

    recorder.processEventsAndGenerateAlerts(
      (3 of uninterestingWriteEvent(workflowId1)) ++
        (5 of uninterestingWriteEvent(workflowId2)) ++
        (7 of uninterestingWriteEvent(workflowId3))
    ) should be(Vector.empty)

    recorder.processEventsAndGenerateAlerts(
      (3 of uninterestingWriteEvent(workflowId1)) ++
        (5 of uninterestingWriteEvent(workflowId2)) ++
        (7 of uninterestingWriteEvent(workflowId3))
    ).toSet should be(Set(HeavyMetadataAlert(workflowId2, 10), HeavyMetadataAlert(workflowId3, 14)))

    recorder.processEventsAndGenerateAlerts(
      (3 of uninterestingWriteEvent(workflowId1)) ++
        (5 of uninterestingWriteEvent(workflowId2)) ++
        (7 of uninterestingWriteEvent(workflowId3))
    ) should be(Vector.empty)

    recorder.processEventsAndGenerateAlerts(
      (3 of uninterestingWriteEvent(workflowId1)) ++
        (5 of uninterestingWriteEvent(workflowId2)) ++
        (7 of uninterestingWriteEvent(workflowId3))
    ).toSet should be(Set(HeavyMetadataAlert(workflowId1, 12), HeavyMetadataAlert(workflowId2, 20), HeavyMetadataAlert(workflowId3, 28)))
  }

  it should "be able to accumulate counts from subworkflows" in {
    val recorder = new ActiveMetadataStatisticsRecorder(10, 10, bundleSubworkflowsIntoParents = true)
    val rootWorkflowId = WorkflowId(UUID.randomUUID())
    val subWorkflow1Id = WorkflowId(UUID.randomUUID())
    val subWorkflow2Id = WorkflowId(UUID.randomUUID())

    // This one is easy to get lost in, hence the liberal comments and scattering of clues

    // Recording SW1 parentage adds 2 events against subWorkflowId1 (and thus rootWorkflowId)
    withClue(recorder.statusString()) {
      recorder.processEventsAndGenerateAlerts(parentNotificationEvent(rootWorkflowId, rootWorkflowId, subWorkflow1Id)) should be(Vector.empty)
    }

    // Recording SW2 parentage adds 2 events against subWorkflowId2 (and thus subWorkflow1Id and thus rootWorkflowId)
    withClue(recorder.statusString()) {
      recorder.processEventsAndGenerateAlerts(parentNotificationEvent(rootWorkflowId, subWorkflow1Id, subWorkflow2Id)) should be(Vector.empty)
    }

    // To get started, add 7 events to the root subworkflow:
    withClue(recorder.statusString()) {
      recorder.processEventsAndGenerateAlerts(7 of uninterestingWriteEvent(rootWorkflowId)) should be(Vector(HeavyMetadataAlert(rootWorkflowId, 11)))
    }

    // Current standing: root: 11, sub1: 4, sub2: 2

    withClue(recorder.statusString()) {
      recorder.processEventsAndGenerateAlerts(7 of uninterestingWriteEvent(subWorkflow1Id)) should be(Vector(HeavyMetadataAlert(subWorkflow1Id, 11)))
    }

    // Current standing: root: 18, sub1: 11, sub2: 2

    withClue(recorder.statusString()) {
      recorder.processEventsAndGenerateAlerts(7 of uninterestingWriteEvent(subWorkflow2Id)) should be(Vector(HeavyMetadataAlert(rootWorkflowId, 25)))
    }

    // Current standing: root: 25, sub1: 18, sub2: 9

    withClue(recorder.statusString()) {
      recorder.processEventsAndGenerateAlerts(7 of uninterestingWriteEvent(subWorkflow1Id)) should be(Vector(HeavyMetadataAlert(subWorkflow1Id, 25)))
    }

    // Current standing: root: 32, sub1: 25, sub2: 9

    withClue(recorder.statusString()) {
      recorder.processEventsAndGenerateAlerts(7 of uninterestingWriteEvent(subWorkflow2Id)).toSet should be(Set(HeavyMetadataAlert(rootWorkflowId, 39), HeavyMetadataAlert(subWorkflow2Id, 16)))
    }

    // Current standing: root: 39, sub1: 32, sub2: 16
  }

  it should "not accumulate counts from subworkflows if disabled" in {
    val recorder = new ActiveMetadataStatisticsRecorder(10, 10, bundleSubworkflowsIntoParents = false)
    val rootWorkflowId = WorkflowId(UUID.randomUUID())
    val subWorkflow1Id = WorkflowId(UUID.randomUUID())
    val subWorkflow2Id = WorkflowId(UUID.randomUUID())

    recorder.processEventsAndGenerateAlerts(parentNotificationEvent(rootWorkflowId, rootWorkflowId, subWorkflow1Id)) should be(Vector.empty)
    recorder.processEventsAndGenerateAlerts(parentNotificationEvent(rootWorkflowId, subWorkflow1Id, subWorkflow2Id)) should be(Vector.empty)

    // If we were accumulating these would alert, but we see nothing if not accumulating:
    recorder.processEventsAndGenerateAlerts(7 of uninterestingWriteEvent(subWorkflow1Id)) should be(Vector.empty)
    recorder.processEventsAndGenerateAlerts(7 of uninterestingWriteEvent(subWorkflow2Id)) should be(Vector.empty)

    // When we trip the limits, we should only see alerts for individual workflows.
    // Note: it's 16 not 14 because of the two parent notification entries above
    recorder.processEventsAndGenerateAlerts(7 of uninterestingWriteEvent(subWorkflow1Id)) should be(Vector(HeavyMetadataAlert(subWorkflow1Id, 16)))
    recorder.processEventsAndGenerateAlerts(7 of uninterestingWriteEvent(subWorkflow2Id)) should be(Vector(HeavyMetadataAlert(subWorkflow2Id, 16)))
  }
}
