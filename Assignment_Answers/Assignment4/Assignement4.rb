require 'bio'
require 'stringio'
require 'fileutils'

#load the documents
protfile= "pepa.fa".to_s
nucfile=  "TAIR10.fa".to_s

#create a directory where they will be located the databases,
#if it exists, delete it and create another and empty one
if File.directory?("./Databases")
  FileUtils.rm_rf "./Databases"
  FileUtils.mkdir_p "./Databases"
end

#create the output file and write the first line explaining what it is,
#if it already exists, delete it and create another empty one
if File.file?("orthologs.txt")
    File.delete("orthologs.txt") #If the file already exists, we delete it
end
orthfile= File.new("orthologs.txt", "a+")
  orthfile.write("This are the orthologue pairs between species Arabidopsis and S. pombe\n\n\n\n") #write the first line with info of the file
orthfile.close

#create both databases
system("makeblastdb -in '#{protfile}' -dbtype 'prot' -out ./Databases/#{protfile}")
system("makeblastdb -in '#{nucfile}' -dbtype 'nucl' -out ./Databases/#{nucfile}")

#create the factories
protFactory= Bio::Blast.local('blastx', "./Databases/#{protfile}")
nucFactory= Bio::Blast.local('tblastn', "./Databases/#{nucfile}")

#and open the files in fasta format to be able to make the blast
protfasta = Bio::FastaFormat.open(protfile)
nucfasta=  Bio::FastaFormat.open(nucfile)

orthlist=[] #to store all the ortholog pairs
cont = 0 #to know in the end how much orthologs are there

#we are going to select the orthologs by evalue and coverage
Evalue = 10**-10 #reference: https://www.biostars.org/p/116338/
cov_threshold= 0.5
fastaseq={}#dict with the ids and secuences 
nucfasta.each do |nucseq|

  fastaseq[(nucseq.entry_id).to_s] = (nucseq.seq).to_s
end

#for each entry, take the queries
protfasta.each {|entry|
  $stderr.puts "\nSearching ... " + entry.definition
  
  protreport = nucFactory.query(entry)
  
  if protreport.hits[0] #if it has any first hit
   
    if protreport.hits[0].evalue <= Evalue #and the evalue lower than the threshold
      
      evalue = protreport.hits[0].evalue #create a variable with the value
      
      #and find the coverage (end-start of the alignment)/the whole sequence
      coverage= (protreport.hits[0].query_end.to_f - protreport.hits[0].query_start.to_f)/protreport.hits[0].query_len.to_f
      
      if  coverage >= cov_threshold #if it is greater than the threshold
        
        #If the coverage and the evalue are above the threshold returns the id of the first hit
        protfind_id = (entry.entry_id).delete("\n").to_s
        protid= protreport.hits[0].definition.split("|")[0].delete(" ").to_s
          
        protq=fastaseq[protid]
 
        #and do the reverse process
        nucreport = protFactory.query(protq)
        
        if nucreport.hits[0]
          #If the coverage and the evalue are above the threshold returns the id of the first hit
          nucprot_id = nucreport.hits[0].definition.split("|")[0].delete(" ").to_s

          if nucprot_id==protfind_id #if both coinccide store in the list and write into the file
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

#the following is for add a message with the total count of ortholog pairs in the third line
orthfile=File.open("orthologs.txt", "a+")
data= orthfile.readlines
orthfile.close

p data

data.insert(2, "There are #{cont} pairs of orthologs")

p data

File.write("orthologs.txt", data.join, mode: "w") # write(overwrite) data to original file