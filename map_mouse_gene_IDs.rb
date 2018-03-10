## ensembl ID map:
# ENSMUSG00000064372      ENSMUST00000082423      mitochondrially encoded tRNA proline [Source:MGI Symbol;Acc:MGI:102478] mt-Tp

## Kallisto output:
# target_id       length  eff_length      est_counts      tpm

def main
  input=ARGV[0]
  mapping = ARGV[1]
  outprefix = ARGV[2]
  if outprefix == nil
    outprefix = input + ".genes"
  end
  transcripts = readMap(mapping)
  
  addGeneNameID(input, transcripts, outprefix)


end

def addGeneNameID(input, transcripts, outprefix)
  genes = {}
  tout = File.new(outprefix + ".transcripts.txt", "w")
  gout = File.new(outprefix + ".genes.txt", "w")

  File.new(input, "r").each do |line|
    if line =~ /^target/  # header
      tout.puts "#{line.chomp}\tGeneID\tGeneName"
    else
      cols = line.chomp.split(/\s+/)
      tid = cols[0].split('.')[0]
      if transcripts.key?(tid)
        gid = transcripts[tid][:geneID]
        gname = transcripts[tid][:geneName]
        if genes.key?(gname) 
          genes[gname] += cols[4].to_f
        else
          genes[gname] = cols[4].to_f
        end
      else
        gid = "NA"
        gname = "NA"
      end
      tout.puts  "#{line.chomp}\t#{gid}\t#{gname}"
    end
  end

  gout.puts "geneName\tTPM"
  genes.keys.sort.each do |gname|
    gout.puts "#{gname}\t#{genes[gname]}"
  end
  tout.close
  gout.close
end 


def readMap(mapping)
  transcripts = {}

  File.new(mapping, "r").each do |line|
    next if line=~ /^\#/  # header lines
    cols = line.chomp.split(/\t/)
    geneID, transcriptID, geneDesc , geneName = cols[0], cols[1], cols[2], cols[3]
    transcripts[transcriptID] = {} unless transcripts.key?(transcriptID)
    
    transcripts[transcriptID][:geneID] = geneID
    transcripts[transcriptID][:geneName] = geneName
    transcripts[transcriptID][:geneDesc] = geneDesc

  end

  return transcripts
end


main()
