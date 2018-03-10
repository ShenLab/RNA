## 
flist = ARGV[0]

outcmdf = flist + ".sh"
outcmd = File.new(outcmdf, "w")

outcmd.puts "#!/bin/bash"

temp = flist + ".concat.txt"
header = "geneName"
i = 1
cmd=""

File.new(flist, "r").each do |f|
#  puts f
  f.chomp!
  sample=f.split('/')[0].sub("kallisto.", "")
  if i > 1
    cmd += " <(cut -f2 #{f})"
  else
    cmd = "paste \"#{f}\""
  end
  header += "\t#{sample}"
  i += 1
end

cmd += "  > #{temp}"
outcmd.puts cmd
## `#{cmd}`

rcmd="sed -i.back \"1 s/^.*/#{header}/\" #{temp}"

outcmd.puts rcmd

outcmd.close

`bash #{outcmdf}`

`rm -f #{temp}.back`
