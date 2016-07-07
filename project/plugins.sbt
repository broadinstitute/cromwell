addSbtPlugin("com.eed3si9n" % "sbt-assembly" % "0.14.1")
addSbtPlugin("io.spray" % "sbt-revolver" % "0.7.2")
addSbtPlugin("org.scalastyle" %% "scalastyle-sbt-plugin" % "0.7.0")
/*
sbt-git 0.7.1 is the last version <= 0.8.5 that works with detached git submodules and our docker build.
See https://github.com/broadinstitute/cromwell/issues/645
 */
addSbtPlugin("com.typesafe.sbt" % "sbt-git" % "0.7.1")
addSbtPlugin("com.github.gseitz" % "sbt-release" % "0.8.3")
addSbtPlugin("org.scoverage" % "sbt-scoverage" % "1.3.5")
addSbtPlugin("org.scoverage" % "sbt-coveralls" % "1.1.0")
