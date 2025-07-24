version 1.0

workflow MemoryRetryTest {
    call TestJavaMemRetry
    call TestJavaKilledRetry
    call TestPythonMemRetry
    call TestOutOfMemoryRetry
    call TestBadCommandRetryMoreMem
}

task TestJavaMemRetry {
    # This is a bit tricky: attempting to allocate an enormous array fails instantly because it *would be* too large, but
    # apparently not because memory was actually exhausted:
    #
    # % java Mem.java
    # Exception in thread "main" java.lang.OutOfMemoryError: Requested array size exceeds VM limit
    #	at Mem.main(Mem.java:1)
    # % echo $?
    # 1
    #
    # But because of the "OutOfMemoryError" that should wind up in stderr, the task would be expected to be retried with
    # more memory, though this may not help.
    command <<<
        echo "MEM_SIZE=$MEM_SIZE" >&2
        echo "MEM_UNIT=$MEM_UNIT" >&2
        echo \
        'class Mem { public static void main(String[] args) { new byte[(int)Math.pow(2, 34)].hashCode(); } }' \
        > Mem.java
        java -Xmx${MEM_SIZE%.*}${MEM_UNIT%?} Mem.java
    >>>
    runtime {
        docker: "eclipse-temurin:21"
        memory: "1 GB"
        maxRetries: 1
    }
}

task TestJavaKilledRetry {
    # This task sets up a -Xmx64g Java VM on a 1 GB memory VM specification plus a memory-hungry workload which should
    # result in it being OOM-killed and retried with more memory.
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

task TestPythonMemRetry {
    # This task allocates a large amount of memory in Python, which should result in its being OOM killed with ax
    # subsequent retry with more memory.
    command <<<
        echo "MEM_SIZE=$MEM_SIZE" >&2
        echo "MEM_UNIT=$MEM_UNIT" >&2
        python3 -c 'print(len([0] * (2**34)))'
    >>>
    runtime {
        docker: "google/cloud-sdk:slim"
        memory: "1 GB"
        maxRetries: 1
    }
}

task TestOutOfMemoryRetry {
    # This task drinks from the /dev/zero firehose, which should result in its being OOM-killed with a subsequent retry
    # with more memory.
    command <<<
        echo "MEM_SIZE=$MEM_SIZE" >&2
        echo "MEM_UNIT=$MEM_UNIT" >&2
        tail /dev/zero
    >>>
    runtime {
        docker: "ubuntu:latest"
        memory: "1 GB"
        maxRetries: 1
    }
}

task TestBadCommandRetryMoreMem {
    # This task runs a command that does not exist on the image. It was also intending to fake its own memory-starvation
    # by writing "Killed" to stderr, but "Killed" has been removed from the retry-with-more-memory configuration. This
    # should simply fail and not be retried.
    command <<<
        echo "MEM_SIZE=$MEM_SIZE" >&2
        echo "MEM_UNIT=$MEM_UNIT" >&2
        echo "Killed" >&2
        bedtools intersect nothing with nothing
    >>>
    runtime {
        docker: "ubuntu:latest"
        memory: "1 GB"
        maxRetries: 1
    }
}

task TestBadCommandRetrySameMem {
    # A task that runs a command that does not exist on the image, and does not attempt to fake its own
    # memory-starvation. This should simply fail and not be retried.
    command <<<
        echo "MEM_SIZE=$MEM_SIZE" >&2
        echo "MEM_UNIT=$MEM_UNIT" >&2
        bedtools intersect nothing with nothing
    >>>
    runtime {
        docker: "ubuntu:latest"
        memory: "1 GB"
        maxRetries: 1
    }
}
