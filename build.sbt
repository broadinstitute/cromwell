import Settings._

resolvers += Resolver.sonatypeRepo("releases")

addCompilerPlugin("org.spire-math" %% "kind-projector" % "0.9.4")

lazy val wom = (project in file("wom"))
  .settings(womSettings: _*)

lazy val wdl = (project in file("wdl"))
  .settings(wdlSettings: _*)
  .dependsOn(wom)

lazy val cwl = (project in file("cwl"))
  .settings(cwlSettings: _*)
  .dependsOn(wom)
  .dependsOn(wom % "test->test")

lazy val root = (project in file("."))
  .settings(rootSettings: _*)
  .aggregate(wom)
  .aggregate(wdl)
  .aggregate(cwl)
  .dependsOn(wdl)
  .dependsOn(cwl)
