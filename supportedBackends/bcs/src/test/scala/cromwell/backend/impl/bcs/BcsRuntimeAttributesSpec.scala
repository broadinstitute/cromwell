package cromwell.backend.impl.bcs

//import cromwell.core.path.DefaultPathBuilder
import wom.values._

class BcsRuntimeAttributesSpec extends BcsTestUtilSpec {
  behavior of s"BcsRuntimeAttributes"

  it should "build correct default runtime attributes from config string" in {
    val runtime = Map.empty[String, WomValue]
    val defaults = createBcsRuntimeAttributes(runtime)
    defaults shouldEqual expectedRuntimeAttributes
  }

  it should "parse dockerTag without docker path" in {
    val runtime = Map("dockerTag" -> WomString("ubuntu/latest"))
    val expected = expectedRuntimeAttributes.copy(dockerTag = Some(BcsDockerWithoutPath("ubuntu/latest")))
    createBcsRuntimeAttributes(runtime) shouldEqual(expected)
  }

  it should "parse dockerTag with path" in {
    val runtime = Map("dockerTag" -> WomString("centos/latest oss://bcs-dir/registry/"))
    val expected = expectedRuntimeAttributes.copy(dockerTag = Some(BcsDockerWithPath("centos/latest", "oss://bcs-dir/registry/")))
    createBcsRuntimeAttributes(runtime) shouldEqual(expected)
  }

  it should "parse dockerTag fail if an empty string value" in {
    val runtime = Map("dockerTag" -> WomString(""))
    an [Exception] should be thrownBy createBcsRuntimeAttributes(runtime)
  }

  it should "parse docker" in {
    val runtime = Map("docker" -> WomString("registry.cn-beijing.aliyuncs.com/test/testubuntu:0.2"))
    val expected = expectedRuntimeAttributes.copy(docker = Some(BcsDockerWithoutPath("registry.cn-beijing.aliyuncs.com/test/testubuntu:0.2")))
    createBcsRuntimeAttributes(runtime) shouldEqual(expected)
  }

  it should "parse correct user data" in {
    val runtime = Map("userData" -> WomString("key value1"))
    val expected = expectedRuntimeAttributes.copy(userData = Some(Vector(BcsUserData("key", "value1"))))
    createBcsRuntimeAttributes(runtime) shouldEqual(expected)
  }

  it should "throw if user data is invalid" in {
    val runtime = Map("userData" -> WomString("keyvalue"))
    an [Exception] should be thrownBy createBcsRuntimeAttributes(runtime)
  }

  it should "parse correct input mount" in {
    val runtime = Map("mounts" -> WomString("oss://bcs-dir/bcs-file /home/inputs/input_file false"))
    //val expected = expectedRuntimeAttributes.copy(mounts = Some(Vector(BcsInputMount(Left(mockPathBuiler.build("oss://bcs-dir/bcs-file").get), Left(DefaultPathBuilder.build("/home/inputs/input_file").get), false))))
    val expected = expectedRuntimeAttributes.copy(mounts = Some(Vector(BcsInputMount(Left(mockPathBuiler.build("oss://bcs-dir/bcs-file").get), Right("/home/inputs/input_file"), false))))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "parse correct out mount" in {
    val runtime = Map("mounts" -> WomString("/home/outputs/ oss://bcs-dir/outputs/ true"))
    //val expected = expectedRuntimeAttributes.copy(mounts = Some(Vector(BcsOutputMount(Left(DefaultPathBuilder.build("/home/outputs/").get), Left(mockPathBuiler.build("oss://bcs-dir/outputs/").get),  true))))
    val expected = expectedRuntimeAttributes.copy(mounts = Some(Vector(BcsOutputMount(Right("/home/outputs/"), Left(mockPathBuiler.build("oss://bcs-dir/outputs/").get),  true))))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "throw if mounts is invalid" in {
    val runtime = Map("mounts" -> WomString("invalid mounts"))
    an [Exception] should be thrownBy createBcsRuntimeAttributes(runtime)
  }

  it should "parse correct cluster id" in {
    val runtime = Map("cluster" -> WomString("cls-1"))
    val expected = expectedRuntimeAttributes.copy(cluster = Some(Left("cls-1")))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "parse correct ondemand auto cluster configuration" in {
    val runtime = Map("cluster" -> WomString("OnDemand ecs.s1.large img-ubuntu"))
    val expected = expectedRuntimeAttributes.copy(cluster = Some(Right(AutoClusterConfiguration("OnDemand", "ecs.s1.large", "img-ubuntu"))))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "parse correct spot auto cluster configuration" in {
    val runtime = Map("cluster" -> WomString("Spot ecs.s1.large img-ubuntu"))
    val expected = expectedRuntimeAttributes.copy(cluster = Some(Right(AutoClusterConfiguration("Spot", "ecs.s1.large", "img-ubuntu"))))
    createBcsRuntimeAttributes(runtime) shouldEqual expected

  }

  it should "parse correct spot auto cluster price option" in {
    val runtime = Map("cluster" -> WomString("Spot ecs.s1.large img-ubuntu SpotWithPriceLimit 0.1"))
    val expected = expectedRuntimeAttributes.copy(cluster = Some(Right(AutoClusterConfiguration("Spot", "ecs.s1.large", "img-ubuntu", Some("SpotWithPriceLimit"), Some(0.1.toFloat)))))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "parse correct vpc cidr block" in {
    val runtime = Map("vpc" -> WomString("172.16.16.0/20"))
    val expected = expectedRuntimeAttributes.copy(vpc = Some(BcsVpcConfiguration(Some("172.16.16.0/20"))))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "parse correct vpc id" in {
    val runtime = Map("vpc" -> WomString("vpc-xxxx"))
    val expected = expectedRuntimeAttributes.copy(vpc = Some(BcsVpcConfiguration(vpcId = Some("vpc-xxxx"))))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "parse correct system disk" in {
    val runtime = Map("systemDisk" -> WomString("cloud_efficiency 250"))
    val expected = expectedRuntimeAttributes.copy(systemDisk = Some(BcsSystemDisk("cloud_efficiency", 250)))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "throw when parsing invalid system disk" in {
    val runtime = Map("systemDisk" -> WomString("cloud_efficiency 250 /home/data/"))
    an [Exception] should be thrownBy createBcsRuntimeAttributes(runtime)
  }

  it should "parse correct data disk" in {
    val runtime = Map("dataDisk" -> WomString("cloud 400 /home/data/"))
    val expected = expectedRuntimeAttributes.copy(dataDisk = Some(BcsDataDisk("cloud", 400, "/home/data/")))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "throw when parsing invalid data disk" in {
    val runtime = Map("dataDisk" -> WomString("cloud_efficiency 250"))
    an [Exception] should be thrownBy createBcsRuntimeAttributes(runtime)
  }

  it should "parse correct reserve on fail option" in {
    val runtime = Map("reserveOnFail" -> WomBoolean(false))
    val expected = expectedRuntimeAttributes.copy(reserveOnFail = Some(false))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "parse correct auto release option" in {
    val runtime = Map("autoReleaseJob" -> WomBoolean(false))
    val expected = expectedRuntimeAttributes.copy(autoReleaseJob = Some(false))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "parse correct verbose option" in {
    val runtime = Map("verbose" -> WomBoolean(false))
    val expected = expectedRuntimeAttributes.copy(verbose = Some(false))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

  it should "parse correct time out" in {
    val runtime = Map("timeout" -> WomInteger(3000))
    val expected = expectedRuntimeAttributes.copy(timeout = Some(3000))
    createBcsRuntimeAttributes(runtime) shouldEqual expected
  }

}
