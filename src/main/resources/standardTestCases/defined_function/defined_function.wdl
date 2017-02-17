task cat {
  File f
  command { cat ${f} }
  runtime { docker: "ubuntu:latest" }
  output { String out = read_string(stdout()) }
}

task mk_file {
  Int index
  command { echo "file_${index}" > out }
  runtime { docker: "ubuntu:latest" }
  output { File f = "out" }
}

# Test the defined function by perhaps calling a function and only following up if the result was
# defined.
workflow defined_function {
  Array[Boolean] mask = [true, false]
  Array[Int] indices = range(2)
  Array[Pair[Boolean, Int]] masked_indices = zip(mask, indices)

  scatter (p in masked_indices) {
    if (p.left) {
      Int good_index = p.right
      call mk_file { input: index = good_index }
    }
    if (defined(mk_file.f)) {
      File real_f = select_first([mk_file.f])
      call cat { input: f = real_f }
    }
  }
  output {
    Array[String?] cat_out = cat.out
    Array[Int] good_indices = select_all(good_index)
  }
}
