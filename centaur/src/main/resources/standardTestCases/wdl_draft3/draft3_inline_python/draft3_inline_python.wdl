version 1.0

workflow draft3_inline_python {
  input {
    Int limit = 100
  }

  call prime_sieve { input: limit = limit }

  output {
    Array[Int] primes = prime_sieve.primes
    Int prime_count = length(primes)
  }
}

task prime_sieve {
  input {
    Int limit
  }

  Int half_limit = floor(limit / 2)

  command <<<
    python <<CODE
    import sys
    import math

    sieve = [True for i in range(~{limit})]
    for i in range(2, ~{half_limit}):
      if sieve[i]:
        for j in range(i * 2, ~{limit}, i):
          sieve[j] = False

    for i in range(2, ~{limit}):
      if sieve[i]:
        print(str(i))

    CODE
  >>>

  runtime {
    docker: "python:latest"
  }

  output {
    Array[Int] primes = read_lines(stdout())
  }
}
