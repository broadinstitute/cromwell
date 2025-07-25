version 1.0

workflow retry_with_more_memory_assorted_ooms {
    call run_app
}

# A Java-based task that sets the maximum heap size to 64 GB on a VM that has nowhere near that much memory. This task
# should be OOM killed and retried with more memory.
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
        java -Xms64g -Xmx64g Mem.java
    >>>
    runtime {
        docker: "eclipse-temurin:21"
        memory: "1 GB"
        maxRetries: 1
    }
}