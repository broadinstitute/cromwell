version 1.0
import "echo_task.wdl" as lib
workflow echo {
    input {
        Array[String] ss
    }
    scatter (s in ss) {
        call lib.echo { input:
            s = s
        }
    }
}
