import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtassembly.{MergeStrategy, PathList}

object Merging {
  val customMergeStrategy: Def.Initialize[String => MergeStrategy] = Def.setting {
    case PathList(ps@_*) if Set("project.properties", "execution.interceptors").contains(ps.last) =>
      // Merge/Filter files from AWS/Google jars that otherwise collide at merge time.
      MergeStrategy.filterDistinctLines
    case PathList(ps@_*) if ps.last == "logback.xml" =>
      MergeStrategy.first
    // AWS SDK v2 configuration files - can be discarded
    case PathList(ps@_*) if Set("codegen.config" , "service-2.json" , "waiters-2.json" , "customization.config" , "examples-1.json" , "paginators-1.json").contains(ps.last) =>
      MergeStrategy.discard
    case x@PathList("META-INF", path@_*) =>
      path map {
        _.toLowerCase
      } match {
        case "spring.tooling" :: _ =>
          MergeStrategy.discard
        case "io.netty.versions.properties" :: Nil =>
          MergeStrategy.first
        case "maven" :: "com.google.guava" :: _ =>
          MergeStrategy.first
        case "versions" :: _ if path.last == "module-info.class" =>
          MergeStrategy.discard
        case "native-image" :: _ if Set("native-image.properties", "reflection-config.json").contains(path.last) =>
          /*
          Discard GraalVM configuration files.
          grpc-netty-shaded 1.39.0 tried to put the netty classes into a different package, but left the shaded version
          of the config file with the same name as the unshaded netty library. Thus when merging the shaded and
          unshaded netty jars we end up with assembly conflicts.

          However, we're not using GraalVM for execution so just discard the configuration files.

          See also:
          - https://www.graalvm.org/reference-manual/native-image/BuildConfiguration/#configuration-file-format
          - https://github.com/grpc/grpc-java/issues/7540
          - https://github.com/grpc/grpc-java/releases/tag/v1.39.0
           */
          MergeStrategy.discard
        case _ =>
          val oldStrategy = (assembly / assemblyMergeStrategy).value
          oldStrategy(x)
      }
    case x@PathList("OSGI-INF", path@_*) =>
      path map {
        _.toLowerCase
      } match {
        case "l10n" :: "bundle.properties" :: Nil =>
          MergeStrategy.concat
        case _ =>
          val oldStrategy = (assembly / assemblyMergeStrategy).value
          oldStrategy(x)
      }
    case "asm-license.txt" | "module-info.class" | "overview.html" | "cobertura.properties" =>
      MergeStrategy.discard
    // inspired by https://github.com/ergoplatform/explorer-backend/blob/7364ecfdeabeb691f0f25525e577d6c48240c672/build.sbt#L14-L15
    case other if other.contains("scala/annotation/nowarn.class")  => MergeStrategy.discard
    case other if other.contains("scala/annotation/nowarn$.class") => MergeStrategy.discard
    case PathList("mime.types") =>
      MergeStrategy.last
    case x =>
      val oldStrategy = (assembly / assemblyMergeStrategy).value
      oldStrategy(x)
  }
}
