package cromwell

import akka.actor.{Props, ActorRef, ActorSystem}
import akka.pattern.ask
import akka.testkit._
import akka.util.Timeout
import better.files.File
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.CromwellTestkitSpec._
import cromwell.backend.{BackendJobDescriptor, BackendConfigurationDescriptor, BackendWorkflowDescriptor, BackendLifecycleActorFactory}
import cromwell.core.WorkflowId
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.backend.{BackendConfigurationEntry, CallLogs}
import cromwell.engine.workflow.OldStyleWorkflowManagerActor
import cromwell.engine.workflow.OldStyleWorkflowManagerActor.{CallStdoutStderr, WorkflowStdoutStderr}
import cromwell.engine.{WorkflowOutputs, _}
import cromwell.server.WorkflowManagerSystem
import cromwell.util.SampleWdl
import cromwell.webservice.CromwellApiHandler._
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, Matchers, OneInstancePerTest, WordSpecLike}
import wdl4s.Call
import wdl4s.values.{WdlArray, WdlFile, WdlString, WdlValue}

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext}
import scala.language.postfixOps
import scala.reflect.ClassTag
import scala.util.matching.Regex

case class TestBackendLifecycleActorFactory(config: Config) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call],
                                                configurationDescriptor: BackendConfigurationDescriptor): Option[Props] = None

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      configurationDescriptor: BackendConfigurationDescriptor): Props = {
    throw new NotImplementedError("this is not implemented")
  }

  override def workflowFinalizationActorProps(): Option[Props] = None
}

object CromwellTestkitSpec {
  val ConfigText =
    """
      |akka {
      |  loggers = ["akka.testkit.TestEventListener"]
      |  loglevel = "INFO"
      |  actor {
      |    debug {
      |       receive = on
      |    }
      |  }
      |  test {
      |    # Some of our tests fire off a message, then expect a particular event message within 3s (the default).
      |    # Especially on CI, the metadata test does not seem to be returning in time. So, overriding the timeouts
      |    # with slightly higher values. Alternatively, could also adjust the akka.test.timefactor only in CI.
      |    filter-leeway = 5s
      |    single-expect-default = 5s
      |    default-timeout = 10s
      |  }
      |}
    """.stripMargin

  val timeoutDuration = 30 seconds

  class TestWorkflowManagerSystem extends WorkflowManagerSystem {
    override protected def systemName: String = "test-system"
    override protected def newActorSystem() = ActorSystem(systemName, ConfigFactory.parseString(CromwellTestkitSpec.ConfigText))
    /**
      * Do NOT shut down the test actor system inside the normal flow.
      * The actor system will be externally shutdown outside the block.
      */
    override def shutdownActorSystem() = {}

    def shutdownTestActorSystem() = super.shutdownActorSystem()
  }

  /**
    * Loans a test actor system. NOTE: This should be run OUTSIDE of a wait block, never within one.
    */
  def withTestWorkflowManagerSystem[T](block: WorkflowManagerSystem => T): T = {
    val testWorkflowManagerSystem = new CromwellTestkitSpec.TestWorkflowManagerSystem
    try {
      block(testWorkflowManagerSystem)
    } finally {
      TestKit.shutdownActorSystem(testWorkflowManagerSystem.actorSystem, timeoutDuration)
    }
  }

  /**
   * Wait for exactly one occurrence of the specified info pattern in the specified block.  The block is in its own
   * parameter list for usage syntax reasons.
   */
  def waitForInfo[T](pattern: String, occurrences: Int = 1)(block: => T)(implicit system: ActorSystem): T = {
    EventFilter.info(pattern = pattern, occurrences = occurrences).intercept {
      block
    }
  }

  /**
   * Wait for occurrence(s) of the specified warning pattern in the specified block.  The block is in its own parameter
   * list for usage syntax reasons.
   */
  def waitForWarning[T](pattern: String, occurrences: Int = 1)(block: => T)(implicit system: ActorSystem): T = {
    EventFilter.warning(pattern = pattern, occurrences = occurrences).intercept {
      block
    }
  }

  /**
   * Wait for occurrence(s) of the specified error pattern in the specified block.  The block is in its own parameter
   * list for usage syntax reasons.
   */
  def waitForError[T](pattern: String, occurrences: Int = 1)(block: => T)(implicit system: ActorSystem): T = {
    EventFilter.error(pattern = pattern, occurrences = occurrences).intercept {
      block
    }
  }

  /**
   * Akka TestKit appears to be unable to match errors generated by `log.error(Throwable, String)` with the normal
   * `EventFilter.error(...).intercept {...}` mechanism since `EventFilter.error` forces the use of a dummy exception
   * that never matches a real exception.  This method works around that problem by building an `ErrorFilter` more
   * explicitly to allow the caller to specify a `Throwable` class.
   */
  def waitForErrorWithException[T](pattern: String,
                                   throwableClass: Class[_ <: Throwable] = classOf[Throwable],
                                   occurrences: Int = 1)
                                  (block: => T)
                                  (implicit system: ActorSystem): T = {
    val regex = Right[String, Regex](pattern.r)
    ErrorFilter(throwableClass, source = None, message = regex, complete = false)(occurrences = occurrences).intercept {
      block
    }
  }

  /**
    * Special case for validating outputs. Used when the test wants to check that an output exists, but doesn't care what
    * the actual value was.
    */
  lazy val AnyValueIsFine: WdlValue = WdlString("Today you are you! That is truer than true! There is no one alive who is you-er than you!")

  implicit class EnhancedWorkflowManagerActor(val manager: ActorRef) extends AnyVal {
    // `implicit` for the asks below.
    private implicit def timeout: Timeout = 30 seconds

    def submit(sources: WorkflowSourceFiles): WorkflowId = {
      val submitMessage = OldStyleWorkflowManagerActor.SubmitWorkflow(sources)
      Await.result(manager.ask(submitMessage), Duration.Inf).asInstanceOf[WorkflowManagerSubmitSuccess].id
    }

    def workflowOutputs(id: WorkflowId): engine.WorkflowOutputs = {
      val outputsMessage = OldStyleWorkflowManagerActor.WorkflowOutputs(id)
      Await.result(manager.ask(outputsMessage).mapTo[WorkflowManagerWorkflowOutputsSuccess], Duration.Inf).outputs
    }

    def workflowStdoutStderr(id: WorkflowId): Map[FullyQualifiedName, Seq[CallLogs]] = {
      val message = WorkflowStdoutStderr(id)
      Await.result(manager.ask(message).mapTo[WorkflowManagerWorkflowStdoutStderrSuccess], Duration.Inf).logs
    }

    def callStdoutStderr(id: WorkflowId, fqn: FullyQualifiedName): Seq[CallLogs] = {
      val message = CallStdoutStderr(id, fqn)
      Await.result(manager.ask(message).mapTo[WorkflowManagerCallStdoutStderrSuccess], Duration.Inf).logs
    }
  }

  lazy val DefaultLocalBackendConfig = ConfigFactory.parseString(
    """
      |  {
      |    root: "cromwell-executions"
      |
      |    filesystems {
      |      local {
      |        localization: [
      |          "hard-link", "soft-link", "copy"
      |        ]
      |      }
      |    }
      |  }
    """.stripMargin)

  lazy val DefaultLocalBackendConfigEntry = BackendConfigurationEntry(
    name = "local",
    lifecycleActorFactoryClass = "cromwell.TestBackendLifecycleActorFactory",
    DefaultLocalBackendConfig
  )

  lazy val JesBackendConfig = ConfigFactory.parseString(
    """
      |{
      |  // Google project
      |  project = "my-cromwell-workflows"
      |
      |  // Base bucket for workflow executions
      |  root = "gs://my-cromwell-workflows-bucket"
      |
      |  // Polling for completion backs-off gradually for slower-running jobs.
      |  // This is the maximum polling interval (in seconds):
      |  maximum-polling-interval = 600
      |
      |  // Optional Dockerhub Credentials. Can be used to access private docker images.
      |  dockerhub {
      |    // account = ""
      |    // token = ""
      |  }
      |
      |  genomics {
      |    // A reference to an auth defined in the `google` stanza at the top.  This auth is used to create
      |    // Pipelines and manipulate auth JSONs.
      |    auth = "service-account"
      |
      |    // Endpoint for APIs, no reason to change this unless directed by Google.
      |    endpoint-url = "https://genomics.googleapis.com/"
      |  }
      |
      |  filesystems {
      |    gcs {
      |      // A reference to a potentially different auth for manipulating files via engine functions.
      |      auth = "user-via-refresh"
      |    }
      |  }
      |}
    """.stripMargin)

  lazy val JesBackendConfigEntry = BackendConfigurationEntry(
    name = "JES",
    lifecycleActorFactoryClass = "cromwell.backend.impl.jes.JesBackendLifecycleActorFactory",
    JesBackendConfig
  )
}

abstract class CromwellTestkitSpec extends TestKit(new CromwellTestkitSpec.TestWorkflowManagerSystem().actorSystem)
  with DefaultTimeout with ImplicitSender with WordSpecLike with Matchers with BeforeAndAfterAll with ScalaFutures with OneInstancePerTest {

  val name = this.getClass.getSimpleName
  implicit val defaultPatience = PatienceConfig(timeout = Span(30, Seconds), interval = Span(100, Millis))
  implicit val ec = system.dispatcher

  def startingCallsFilter[T](callNames: String*)(block: => T): T =
    waitForPattern(s"starting calls: ${callNames.mkString(", ")}$$") {
      block
    }

  def waitForHandledMessage[T](named: String)(block: => T): T = {
    waitForHandledMessagePattern(s"^received handled message $named") {
      block
    }
  }

  def waitForHandledMessagePattern[T](pattern: String)(block: => T): T = {
    waitForInfo(pattern = pattern, occurrences = 1) {
      block
    }
  }

  /**
   * Performs the following steps:
   *
   * <ol>
   * <li> Sends the specified message to the implicitly passed `ActorRef` via an `ask`.
   * <li> Collects the `Future[Any]` response.
   * <li> Downcasts the `Future[Any]` to a `Future[M]`.
   * <li> Issues a blocking `Await.result` on the `Future`, yielding an `M`.
   * </ol>
   *
   */
  def messageAndWait[M: ClassTag](message: AnyRef)(implicit actorRef: ActorRef): M = {
    val futureAny = actorRef ? message
    Await.result(futureAny.mapTo[M], timeoutDuration)
  }

  /**
   * Wait for exactly one occurrence of the specified info pattern in the specified block within a limited amount of
   * time. The block is in its own parameter list for usage syntax reasons.
   */
  def waitForPattern[T](pattern: String, occurrences: Int = 1)(block: => T): T = {
    within(timeoutDuration) {
      waitForInfo(pattern, occurrences) {
        block
      }
    }
  }

  private def buildWorkflowManagerActor(config: Config) = {
    TestActorRef(new OldStyleWorkflowManagerActor(config))
  }

  // Not great, but this is so we can test matching data structures that have WdlFiles in them more easily
  private def validateOutput(output: WdlValue, expectedOutput: WdlValue): Unit = expectedOutput match {
    case expectedFile: WdlFile if output.isInstanceOf[WdlFile] =>
      val actualFile = output.asInstanceOf[WdlFile]
      actualFile.value.toString.endsWith(expectedFile.value.toString) shouldEqual true
    case expectedArray: WdlArray if output.isInstanceOf[WdlArray] =>
      val actualArray = output.asInstanceOf[WdlArray]
      actualArray.value.size should be(expectedArray.value.size)
      (actualArray.value zip expectedArray.value) foreach {
        case (actual, expected) => validateOutput(actual, expected)
      }
    case _ =>
      output shouldEqual expectedOutput
  }

  def runWdl(sampleWdl: SampleWdl,
             eventFilter: EventFilter,
             runtime: String = "",
             workflowOptions: String = "{}",
             terminalState: WorkflowState = WorkflowSucceeded,
             config: Config = OldStyleWorkflowManagerActor.defaultConfig)(implicit ec: ExecutionContext): WorkflowOutputs = {
    val wma = buildWorkflowManagerActor(config)
    val sources = WorkflowSourceFiles(sampleWdl.wdlSource(runtime), sampleWdl.wdlJson, workflowOptions)
    eventFilter.intercept {
      within(timeoutDuration) {
        val workflowId = wma.submit(sources)
        verifyWorkflowState(wma, workflowId, terminalState)
        wma.workflowOutputs(workflowId)
      }
    }
  }

  def runWdlAndAssertOutputs(sampleWdl: SampleWdl,
                             eventFilter: EventFilter,
                             expectedOutputs: Map[FullyQualifiedName, WdlValue],
                             runtime: String = "",
                             workflowOptions: String = "{}",
                             allowOtherOutputs: Boolean = true,
                             terminalState: WorkflowState = WorkflowSucceeded,
                             config: Config = OldStyleWorkflowManagerActor.defaultConfig,
                             workflowManagerActor: Option[TestActorRef[OldStyleWorkflowManagerActor]] = None)
                            (implicit ec: ExecutionContext): WorkflowId = {
    val workflowManager = workflowManagerActor.getOrElse(buildWorkflowManagerActor(config))
    val sources = sampleWdl.asWorkflowSources(runtime, workflowOptions)
    eventFilter.intercept {
      within(timeoutDuration) {
        val workflowId = workflowManager.submit(sources)
        verifyWorkflowState(workflowManager, workflowId, terminalState)
        val outputs = workflowManager.workflowOutputs(workflowId)

        val actualOutputNames = outputs.keys mkString ", "
        val expectedOutputNames = expectedOutputs.keys mkString " "

        expectedOutputs foreach { case (outputFqn, expectedValue) =>
          val actualValue = outputs.getOrElse(outputFqn, throw new RuntimeException(s"Expected output $outputFqn was not found in: '$actualOutputNames'"))
          if (expectedValue != AnyValueIsFine) validateOutput(actualValue.wdlValue, expectedValue)
        }
        if (!allowOtherOutputs) {
          outputs foreach { case (actualFqn, actualValue) =>
            val expectedValue = expectedOutputs.getOrElse(actualFqn, throw new RuntimeException(s"Actual output $actualFqn was not wanted in '$expectedOutputNames'"))
            if (expectedValue != AnyValueIsFine) validateOutput(actualValue.wdlValue, expectedValue)
          }
        }
        workflowId
      }
    }
  }

  /*
     FIXME: I renamed this as it appears to be asserting the stdout/stderr of a single call which is kinda weird for
     a full workflow type of thing
  */
  def runSingleCallWdlWithWorkflowManagerActor(wma: TestActorRef[OldStyleWorkflowManagerActor],
                                               sources: WorkflowSourceFiles,
                                               eventFilter: EventFilter,
                                               fqn: FullyQualifiedName,
                                               index: ExecutionIndex,
                                               stdout: Option[Seq[String]],
                                               stderr: Option[Seq[String]],
                                               expectedOutputs: Map[FullyQualifiedName, WdlValue] = Map.empty)(implicit ec: ExecutionContext) = {
    eventFilter.intercept {
      within(timeoutDuration) {
        val workflowId = wma.submit(sources)
        verifyWorkflowState(wma, workflowId, WorkflowSucceeded)
        val standardStreams = wma.callStdoutStderr(workflowId, fqn)
        stdout foreach { souts =>
          souts shouldEqual (standardStreams map { s => File(s.stdout.value).contentAsString })
        }
        stderr foreach { serrs =>
          serrs shouldEqual (standardStreams map { s => File(s.stderr.value).contentAsString })
        }
      }
    }
  }

  def runWdlWithWorkflowManagerActor(wma: TestActorRef[OldStyleWorkflowManagerActor],
                                     sources: WorkflowSourceFiles,
                                     eventFilter: EventFilter,
                                     stdout: Map[FullyQualifiedName, Seq[String]],
                                     stderr: Map[FullyQualifiedName, Seq[String]],
                                     expectedOutputs: Map[FullyQualifiedName, WdlValue] = Map.empty,
                                     terminalState: WorkflowState = WorkflowSucceeded,
                                     assertStdoutStderr: Boolean = false)(implicit ec: ExecutionContext) = {
    eventFilter.intercept {
      within(timeoutDuration) {
        val workflowId = wma.submit(sources)
        verifyWorkflowState(wma, workflowId, terminalState)

        if (assertStdoutStderr) {
          val standardStreams = wma.workflowStdoutStderr(workflowId)
          stdout foreach {
            case (fqn, out) if standardStreams.contains(fqn) =>
              out shouldEqual (standardStreams(fqn) map { s => File(s.stdout.value).contentAsString })
          }
          stderr foreach {
            case (fqn, err) if standardStreams.contains(fqn) =>
              err shouldEqual (standardStreams(fqn) map { s => File(s.stderr.value).contentAsString })
          }
        }
      }
    }
  }

  def runWdlAndAssertStdoutStderr(sampleWdl: SampleWdl,
                                  eventFilter: EventFilter,
                                  fqn: FullyQualifiedName,
                                  index: ExecutionIndex,
                                  runtime: String = "",
                                  stdout: Option[Seq[String]] = None,
                                  stderr: Option[Seq[String]] = None,
                                  config: Config = OldStyleWorkflowManagerActor.defaultConfig)
                                 (implicit ec: ExecutionContext) = {
    val actor = buildWorkflowManagerActor(config)
    val sources = WorkflowSourceFiles(sampleWdl.wdlSource(runtime), sampleWdl.wdlJson, "{}")
    runSingleCallWdlWithWorkflowManagerActor(actor, sources, eventFilter, fqn, index, stdout, stderr)
  }

  def runWdlAndAssertWorkflowStdoutStderr(sampleWdl: SampleWdl,
                                          eventFilter: EventFilter,
                                          runtime: String = "",
                                          stdout: Map[FullyQualifiedName, Seq[String]] = Map.empty[FullyQualifiedName, Seq[String]],
                                          stderr: Map[FullyQualifiedName, Seq[String]] = Map.empty[FullyQualifiedName, Seq[String]],
                                          terminalState: WorkflowState = WorkflowSucceeded,
                                          config: Config = OldStyleWorkflowManagerActor.defaultConfig)
                                         (implicit ec: ExecutionContext) = {
    val actor = buildWorkflowManagerActor(config)
    val workflowSources = WorkflowSourceFiles(sampleWdl.wdlSource(runtime), sampleWdl.wdlJson, "{}")
    runWdlWithWorkflowManagerActor(actor, workflowSources, eventFilter, stdout, stderr, Map.empty, terminalState)
  }

  def verifyWorkflowState(wma: ActorRef, workflowId: WorkflowId, expectedState: WorkflowState)(implicit ec: ExecutionContext): Unit = {
    // Continuously check the state of the workflow until it is in a terminal state
    awaitCond(pollWorkflowState(wma, workflowId).isTerminal)
    // Now that it's complete verify that we ended up in the state we're expecting
    pollWorkflowState(wma, workflowId) should equal (expectedState)
  }

  private def pollWorkflowState(wma: ActorRef, workflowId: WorkflowId)(implicit ec: ExecutionContext): WorkflowState = {
    val futureResponse = wma.ask(OldStyleWorkflowManagerActor.WorkflowStatus(workflowId)).mapTo[WorkflowManagerStatusSuccess]
    Await.result(futureResponse map { _.state }, timeoutDuration)
  }
}
