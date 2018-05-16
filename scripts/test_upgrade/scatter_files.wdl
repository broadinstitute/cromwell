task mkFile {

 Int index

 command {
   echo "content-${index}"
 }

 output {
  File f = stdout()
 }
 runtime { docker: "ubuntu:latest" }
}

task catFile {
 File f_in

 command {
   sleep 50
   cat ${f_in}
 }

 output {
  File f = stdout()
 }
}

task gather {
 Array[File] inputs

 command {
   cat ${sep=" " inputs}
 }

 output {
  String result = stdout()
 }
 runtime { docker: "ubuntu:latest" }
}

workflow scatter_files {

 Int scatter_width_part_1
 Int scatter_width_part_2
 Int scatter_width_part_3
 Int scatter_width_part_4
 Int scatter_width_part_5
 Int scatter_width_part_6

 Int sum = scatter_width_part_1 + scatter_width_part_2 + scatter_width_part_3 + scatter_width_part_4 + scatter_width_part_5 + scatter_width_part_6

 Array[Int] xs = range(sum / 4)

 scatter (x in xs) {
  call mkFile { input: index = x }
  call catFile { input: f_in = mkFile.f}
 }

 call gather { input: inputs = catFile.f }
}
