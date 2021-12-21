require 'bio'
require 'stringio'
require 'fileutils'


protfile= "pepa.fa".to_s
nucfile=  "TAIR10.fa".to_s

if File.directory?("./Databases")
  FileUtils.rm_rf "./Databases"
  FileUtils.mkdir_p "./Databases"
end

if File.file?("orthologs.txt")
    File.delete("orthologs.txt") #If the file already exists, we delete it
end
orthfile= File.new("orthologs.txt", "a+")
  orthfile.write("This are the orthologue pairs between species Arabidopsis and S. pombe\n\n\n\n") #write the first line with info of the file
orthfile.close


system("makeblastdb -in '#{protfile}' -dbtype 'prot' -out ./Databases/#{protfile}")
system("makeblastdb -in '#{nucfile}' -dbtype 'nucl' -out ./Databases/#{nucfile}")

protFactory= Bio::Blast.local('blastx', "./Databases/#{protfile}")
nucFactory= Bio::Blast.local('tblastn', "./Databases/#{nucfile}")

protfasta = Bio::FastaFormat.open(protfile)
nucfasta=  Bio::FastaFormat.open(nucfile)

orthlist=[]
cont = 0 #to know in the end how much orthologs are there
Evalue = 10**-10
cov_threshold= 0.5
fastaseq={}
nucfasta.each do |nucseq|

  fastaseq[(nucseq.entry_id).to_s] = (nucseq.seq).to_s
end

protfasta.each {|entry|
  $stderr.puts "\nSearching ... " + entry.definition
  
  protreport = nucFactory.query(entry)
  
  if protreport.hits[0]
   
    if protreport.hits[0].evalue <= Evalue
      
      evalue = protreport.hits[0].evalue
      
      coverage= (protreport.hits[0].query_end.to_f - protreport.hits[0].query_start.to_f)/protreport.hits[0].query_len.to_f
      
      if  coverage >= cov_threshold
      
        protfind_id = (entry.entry_id).delete("\n").to_s
        protid= protreport.hits[0].definition.split("|")[0].delete(" ").to_s
          
        protq=fastaseq[protid]
            
        nucreport = protFactory.query(protq)
        
        if nucreport.hits[0]
          nucprot_id = nucreport.hits[0].definition.split("|")[0].delete(" ").to_s
          
          if nucprot_id==protfind_id
            unless orthlist.include?("#{nucprot_id}\t#{protfind_id}")
              orthlist << "#{nucprot_id}, #{protfind_id}"
              cont +=1
              orthfile=File.open("orthologs.txt", "a+")
                orthfile.write("#{nucprot_id}\t#{protfind_id}, coverage = #{coverage}, evalue = #{evalue}\n")
              orthfile.close

            end

          end
         
        end
        
      end
    
    end
  
  end
    }

orthfile=File.open("orthologs.txt", "a+")
data= orthfile.readlines
orthfile.close

p data

data.insert(2, "There are #{cont} pairs of orthologs")

p data

File.write("orthologs.txt", data.join, mode: "w") # write(overwrite) data to original file