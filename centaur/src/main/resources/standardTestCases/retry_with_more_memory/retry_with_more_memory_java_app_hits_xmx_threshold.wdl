version 1.0

workflow retry_with_more_memory_java_app_hits_xmx_threshold {
    call run_app
}

# A Java-based task that sets the maximum heap size to 1 GB but allocates far more memory than that. This should cause
# an OOM to be written to stderr which combined with a `memory_retry_multiplier` of 2, should produce a second attempt
# of the job with double the memory of the first attempt.
task run_app {
    command <<<
        echo "MEM_SIZE=$MEM_SIZE" >&2
        echo "MEM_UNIT=$MEM_UNIT" >&2

        cat > Mem.java << EOF
        class Mem {
          public static void main(String[] args) throws Exception {
            int gb = (int)Math.pow(2, 30);
            System.out.println("Allocating memory...");
            byte[][] byteArr = new byte[32][];
            for (int i = 0; i < byteArr.length; i++) {
                byteArr[i] = new byte[gb];
            }
            System.out.println("Sleeping a minute...");
            Thread.sleep(60_000);
            System.out.printf("Heap size: %,.2f%n", (double)Runtime.getRuntime().totalMemory() / gb);
            System.out.println(byteArr.hashCode());
          }
        }

        EOF
        java -Xmx1G Mem.java
    >>>
    runtime {
        docker: "eclipse-temurin:21"
        memory: "4 GB"
        maxRetries: 1
    }
}