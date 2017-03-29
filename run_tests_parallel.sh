#
# This approach allows for the running of scalatests that mixin the
# trait 'ParallelTestExecution' to be executed in parallel with a
# throttle to the number of threads.  Simply running 'sbt test' 
# has that throttle hardcoded to be 2*cores
#

# optional parameter for parallelism, defaults to 3
THREADS=${1:-3}

echo "Running tests with ${THREADS}-way parallelism"

sbt test:compile
CP=$(sbt "export test:dependency-classpath" --error)
java -cp $CP org.scalatest.tools.Runner -R target/scala-2.12/test-classes -oD -PS${THREADS}
