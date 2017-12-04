task A {
  command {
    echo -n -e "jeff\nchris\nmiguel\nthibault\nkhalid\nscott"
  }
  output {
    Array[String] A_out = read_lines(stdout())
  }
}

task B {
  String B_in
  command {
    python -c "print(len('${B_in}'))"
  }
  output {
    Int B_out = read_int(stdout())
  }
}

task C {
  Int C_in
  command {
    python -c "print(${C_in}*100)"
  }
  output {
    Int C_out = read_int(stdout())
  }
}

task D {
  Array[Int] D_in
  command {
    python -c "print(${sep = '+' D_in})"
  }
  output {
    Int D_out = read_int(stdout())
  }
}

task E {
  command {
    python -c "import random; print(random.randint(1,100))"
  }
  output {
    Int E_out = read_int(stdout())
  }
}

workflow w {
  call A
  scatter (item in A.A_out) { # scatter 0
    call B {input: B_in = item}
    call C {input: C_in = B.B_out}
    call E
    scatter (itemB in B.B_out) { # scatter 1
      call E as G
    }
    scatter (itemB in B.B_out) { # scatter 2
      call E as H
    }
  }
  scatter (item in A.A_out) { # scatter 3
    call E as F
  }
  call D {input: D_in = B.B_out}
}