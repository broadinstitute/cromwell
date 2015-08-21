package cromwell.engine

import cromwell.binding._
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow._
import cromwell.parser.BackendType
import cromwell.util.SampleWdl
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers}

class ScatterPopulateSpec extends FlatSpec with Matchers {

  import ScatterPopulateSpec._

  behavior of "ScatterPopulate"

  it should "populate WDL from our conversation" in {
    val namespace = NamespaceWithWorkflow.load(SampleWdl.ScatterWdl.wdlSource(), BackendType.LOCAL)

    val workflowEntries = WorkflowActor.populate(namespace.workflow)
    workflowEntries.size should be(3)

    val workflowEntriesExpected = Table(
      ("fqn", "type"),
      ("w.A", classOf[CallKey]),
      ("w.$scatter_0", classOf[ScatterKey]),
      ("w.D", classOf[CallKey]))

    forAll(workflowEntriesExpected) {
      (fqn: FullyQualifiedName, clazz: Class[_]) =>
        validateWorkflow(workflowEntries, fqn, clazz)
    }

    val scatter0 = workflowEntries.findScatterKey("w.$scatter_0")
    val scatter0Entries = scatter0.populate(3)
    scatter0Entries.size should be(12)

    Seq("w.B", "w.C", "w.E") foreach { fqn =>
      validateScatter(scatter0Entries, fqn, classOf[CallKey], scatter0, 3)
    }
  }

  it should "populate scatter WDL example" in {
    val namespace = NamespaceWithWorkflow.load(SampleWdl.NestedScatterWdl.wdlSource(), BackendType.LOCAL)

    // Workflow: Populate the workflow entries
    val workflowEntries = WorkflowActor.populate(namespace.workflow)
    workflowEntries.size should be(4)

    val workflowEntriesExpected = Table(
      ("fqn", "type"),
      ("w.A", classOf[CallKey]),
      ("w.$scatter_0", classOf[ScatterKey]),
      ("w.$scatter_3", classOf[ScatterKey]),
      ("w.D", classOf[CallKey]))

    forAll(workflowEntriesExpected) {
      (fqn: FullyQualifiedName, clazz: Class[_]) =>
        validateWorkflow(workflowEntries, fqn, clazz)
    }

    // From the workflow entries, populate the first w.$scatter_0

    val scatter0 = workflowEntries.findScatterKey("w.$scatter_0")
    val scatter0Entries = scatter0.populate(3)
    scatter0Entries.size should be(20)

    val scatter0EntriesExpected = Table(
      ("fqn", "type"),
      ("w.B", classOf[CallKey]),
      ("w.C", classOf[CallKey]),
      ("w.E", classOf[CallKey]),
      ("w.$scatter_1", classOf[ScatterKey]),
      ("w.$scatter_2", classOf[ScatterKey]))

    forAll(scatter0EntriesExpected) {
      (fqn: FullyQualifiedName, clazz: Class[_]) =>
        validateScatter(scatter0Entries, fqn, clazz, scatter0, 3)
    }

    // From Scatter 0 entries, find Scatter 2.
    // Scatter 0 created three instances of all children, including Scatter 2.
    // Get the child at Index 1, and populate four entries

    val scatter2 = scatter0Entries.findScatterKey("w.$scatter_2", Option(1))
    val scatter2Entries = scatter2.populate(4)
    scatter2Entries.size should be(5)
    validateScatter(scatter2Entries, "w.H", classOf[CallKey], scatter2, 4)
  }

  def validateWorkflow(entries: ExecutionStore, fqn: FullyQualifiedName, clazz: Class[_]): Unit = {
    val key = entries.findKey(fqn, None)
    key.getClass should be(clazz)
    key.parent should be(empty)
    entries.map(_.status) should contain only ExecutionStatus.NotStarted
  }

  def validateScatter(entries: ExecutionStore, fqn: FullyQualifiedName, clazz: Class[_],
                      parent: ScatterKey, count: Int): Unit = {
    val collectorKey = entries.findKey(fqn, None)
    collectorKey shouldBe a[CollectorKey]
    collectorKey.parent shouldNot be(empty)
    collectorKey.parent.get should be(parent)

    (0 until count) foreach { i =>
      val index = Option(i)
      val key = entries.findKey(fqn, index)
      key.getClass should be(clazz)
      key.parent shouldNot be(empty)
      key.parent.get should be(parent)
    }
    entries.map(_.status) should contain only ExecutionStatus.NotStarted
  }
}

object ScatterPopulateSpec {

  implicit class ExecutionStoreWrapper(val entries: ExecutionStore) extends AnyVal {
    def findKeyOption(fqn: FullyQualifiedName, index: Option[Int] = None): Option[ExecutionStoreKey] = {
      entries map {
        _.key
      } find { key =>
        key.scope.fullyQualifiedName == fqn && key.index == index
      }
    }

    def findKey(fqn: FullyQualifiedName, index: Option[Int] = None) = {
      val findOption = findKeyOption(fqn, index)
      if (findOption.isDefined) {
        findOption.get
      } else {
        val keys = entries map { entry =>
          s"${entry.key.scope.fullyQualifiedName} || ${entry.key.index}"
        }
        val msg = keys.mkString(s"$fqn || $index not found in:\n", "\n", "")
        throw new NoSuchElementException(msg)
      }
    }

    def findScatterKey(fqn: FullyQualifiedName, index: Option[Int] = None): ScatterKey = {
      findKey(fqn, index).asInstanceOf[ScatterKey]
    }
  }

  implicit class ExecutionStoreEntryWrapper(val entry: ExecutionStoreEntry) extends AnyVal {
    def key = entry._1

    def status = entry._2

    override def toString = s"ExecutionStoreEntry(key=$key,status=$status)"
  }

  implicit class ExecutionStoreKeyWrapper(val key: ExecutionStoreKey) extends AnyVal {
    def scope = key.scope

    def index = key.index

    def parent = key.parent

    def typeName = key.getClass.getSimpleName.stripSuffix("Key")

    override def toString = s"ExecutionStoreKey(scope=$scope,index=$index,parent=$parent)"
  }

}
