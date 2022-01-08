package cromwell.backend.impl.bcs.callcaching


import akka.actor.Props
import akka.testkit.TestActorRef
import com.typesafe.config.ConfigValueFactory
import cromwell.backend.impl.bcs.{BcsBackendInitializationData, BcsConfiguration, BcsRuntimeAttributes, BcsTestUtilSpec, BcsWorkflowPaths}
import cromwell.backend.standard.callcaching.StandardCacheHitCopyingActorParams
import cromwell.core.path.{Path}
import wom.values._
import cromwell.backend.impl.bcs.BcsTestUtilSpec.BcsBackendConfig
import cromwell.backend.standard.callcaching.DefaultStandardCacheHitCopyingActorParams
import cromwell.core.simpleton.WomValueSimpleton
import cromwell.filesystems.oss.OssPath
import org.mockito.Mockito.when

import scala.util.Try


class BcsBackendCacheHitCopyingActorSpec extends BcsTestUtilSpec  {
  behavior of "BcsBackendCacheHitCopyingActor"
  type ValueOrDelete = Either[Boolean, AnyRef]

  val workflowPaths = BcsWorkflowPaths(workflowDescriptor, BcsBackendConfig, mockPathBuilders)


  private def buildInitializationData(configuration: BcsConfiguration) = {

    val runtimeAttributesBuilder = BcsRuntimeAttributes.runtimeAttributesBuilder(BcsTestUtilSpec.BcsBackendConfigurationDescriptor.backendRuntimeAttributesConfig)
    BcsBackendInitializationData(workflowPaths, runtimeAttributesBuilder, configuration, null)
  }

  val runtimeAttributesBuilder = BcsRuntimeAttributes.runtimeAttributesBuilder(BcsTestUtilSpec.BcsBackendConfigurationDescriptor.backendRuntimeAttributesConfig)


  private def withConfig(configs: Map[String, ValueOrDelete]) = {
    var descriptor = BcsTestUtilSpec.BcsBackendConfigurationDescriptor.copy()
    for ((key, value) <- configs) {
      value match {
        case Left(_) => descriptor = BcsTestUtilSpec.BcsBackendConfigurationDescriptor.copy(backendConfig = descriptor.backendConfig.withoutPath(key))
        case Right(v) => descriptor = BcsTestUtilSpec.BcsBackendConfigurationDescriptor.copy(backendConfig = descriptor.backendConfig.withValue(key, ConfigValueFactory.fromAnyRef(v)))
      }
    }
    new BcsConfiguration(descriptor)
  }


  var cacheHitCopyingActorParams = {
    val mockCacheHitCopyingActorParams = mock[DefaultStandardCacheHitCopyingActorParams]
    val id = "test-access-id"
    val key = "test-access-key"
    val configs = Map("access-id" -> Right(id), "access-key" -> Right(key))
    val conf = withConfig(configs)
    when(mockCacheHitCopyingActorParams.backendInitializationDataOption).thenReturn(Option(buildInitializationData(conf)))
    mockCacheHitCopyingActorParams
  }

  class TestableBcsCacheHitCopyingActor(params: StandardCacheHitCopyingActorParams)
    extends BcsBackendCacheHitCopyingActor(params) {

    val id = "test-access-id"
    val key = "test-access-key"
    val configs = Map("access-id" -> Right(id), "access-key" -> Right(key))
    val conf = withConfig(configs)


    def this() = {
      this(cacheHitCopyingActorParams)
    }

    override def getPath(str: String): Try[Path] = mockPathBuilder.build("oss://bcs-dir/outputs/")
    override def destinationJobDetritusPaths: Map[String, Path] = Map("stdout"
      -> mockPathBuilder.build("oss://my-bucket/cromwell_dir/wf_echo/14e5dcd2-0c94-4035-aa7b-b90d7008202c/call-echo/stdout.log").get)
  }

  it should "process simpleton and detritus correctly" in {
    val simpleton = WomValueSimpleton("txt_files", WomSingleFile("oss://my-bucket/cromwell_dir/wf_echo/14e5dcd2-0c94-4035-aa7b-b90d7008202c/call-echo/abc.log"))
    val detritus = Map("stdout" -> "oss://my-bucket/cromwell_dir/wf_echo/14e5dcd2-0c94-4035-aa7b-b90d7008202c/call-echo/stdout.log")
    val sourceCallRootPath: OssPath = mockPathBuilder.build("oss://bcs-test/root/abc.log").getOrElse(throw new IllegalArgumentException())

    val props = Props(new TestableBcsCacheHitCopyingActor())
    val cacheHitActor = TestActorRef[TestableBcsCacheHitCopyingActor](
      props, "TestableBcsCacheHitCopyingActor")

    noException should be thrownBy cacheHitActor.underlyingActor.processSimpletons(List(simpleton), sourceCallRootPath)
    noException should be thrownBy cacheHitActor.underlyingActor.processDetritus(detritus)
  }
}
