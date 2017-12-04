import "subtest.wdl" as sub

task t {
 String s
 command {...}
 output {
   String o = "yahh"
   Int o2 = 5
 }
}

workflow w {
 call t { input: s = "he" }
 call sub.t as subt { input: s = "ho" }
 call sub.w2 as sw2 { input: s = "ha" }
 call t as u { input: s = sw2.result }
 call sub.w2 as sub2 { input: inp = u.o }
}