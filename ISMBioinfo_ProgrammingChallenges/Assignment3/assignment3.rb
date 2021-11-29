require 'bio'
require 'net/http'
require 'rest-client'

def fetch(url, headers = {accept: "*/*"}, user = "", pass="") #to handle errors
    response = RestClient::Request.execute({
      method: :get,
      url: url.to_s,
      user: user,
      password: pass,
      headers: headers})
    return response
    rescue RestClient::ExceptionWithResponse => e 
      $stderr.puts e.inspect
      response = false
      return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
    rescue RestClient::Exception => e
      $stderr.puts e.inspect
      response = false
      return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
    rescue Exception => e
      $stderr.puts e.inspect
      response = false
      return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  end


def create_features(coord,seq, str, name, motif, attr) # the function needs: coord -> an array with the coordinates of the genes/chromosomes
  #depending on what we are working with;
  #seq-> the biosequence of each gene
  #str -> an array with the sign of the strands where the exons are
  #name -> name of the feature
  #motif -> in this case, cttctt
  #attr -> an array with the attributes (the nomber of the exon)
  #With that, add a new feature to the sequence object 
    coord.each{|c| 
    f1 = Bio::Feature.new(name, c ) 
    f1.append(Bio::Feature::Qualifier.new('repeat_motif', motif))
    f1.append(Bio::Feature::Qualifier.new('note', 'found by repeatfinder'))
    f1.append(Bio::Feature::Qualifier.new('strand', str[coord.index(c)])) 
    f1.append(Bio::Feature::Qualifier.new('attribute', attr[coord.index(c)]))
    f1.append(Bio::Feature::Qualifier.new('from', c.split("..")[0])) #take the first number of the coordinate
    f1.append(Bio::Feature::Qualifier.new('to', c.split("..")[1])) #take the first number of the coordinate
    seq.features << f1  # Append created feature
    } 
  
end

def write_gff(file,geneid,bioseq) # Write a new gff file with the gene where the repeated motif is, given the file name,
  #the genename/chromosome number and a bioseq feature 
  gff = File.open(file, "a+")
  source = "BioRuby"
  type = "cttctt_repeat"
  score = "."
  phase = "."
  gene_id = geneid

    strand = bioseq.assoc["strand"]
    attributes = bioseq.assoc["attribute"]
    from = bioseq.assoc["from"]
    to = bioseq.assoc["to"]
    
    gff.write("#{gene_id}\t#{source}\t#{type}\t#{from}\t#{to}\t#{score}\t#{strand}\t#{phase}\t#{attributes}\n")
end

listfile = File.open("/home/osboxes/Bioinfo_ProgrammingChallenges/Assignment_Answers/Assignment3/ArabidopsisSubNetwork_GeneList.txt")
  if File.file?("gff_gene.gff")
    File.delete("gff_gene.gff") #If the file already exists, we delete it
  end
  gff= File.new("gff_gene.gff", "a+")
    gff.write("##gff-version 3\n\n") #write the first line with info of the file
  gff.close
  
  if File.file?("genes_without_cttctt.txt")
    File.delete("genes_without_cttctt.txt") #If the file already exists, we delete it
  end
  report1= File.new("genes_without_cttctt.txt", "a+")
    report1.write("##This is the list of genes that don't have CTTCTT in their sequence.\n\n")
  report1.close
  
  if File.file?("gff_chr.gff")
    File.delete("gff_chr.gff") #If the file already exists, we delete it
  end
  gff= File.new("gff_chr.gff", "a+")
    gff.write("##gff-version 3\n\n")
  gff.close
  
  arabidopsisgenes = [] #array with all the genes
  
    File.readlines(listfile).each do |line| #read each line of the document
      line=line.strip #remove the \n at the end of each line
      arabidopsisgenes << line
      address = fetch("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{line}")  
      if address #if the adress is valid
                                             
        record = address.body #take the body (data)
        entry = Bio::EMBL.new(record) #create new object
        emblseq = entry.to_biosequence #get the biosequence
        seqembl= entry.seq
        exonlist=[] #array with all the exon sequences
        genecoordlist=[] #array with the coordinates of the genes
        exoncoordlist=[] #array with the coordinates of the exons in the genes
        chrcoordlist= [] #array with the coordinates of the chromosomes
        geneinlist=[] #array with the gene id's that has 'cttctt' motif
        strandlist = [] #array with the signs of the strand where the 'cttctt' motif is
        attributelist=[] #array with the information of the exons that have 'cttctt' motif
        nocttgenes= [] #array with the genes that don't have any 'cttctt' motif
        
        m=entry.accession.match(/chromosome:[A-Z0-9]+:(?<chromosome>\d+):(?<start_chr>\d+):(?<end_chr>\d+)/) #Get chromosome coordinates of gene
        chr= m[:chromosome] #get the number of the chrmosome
        chr_from= m[:start_chr].to_i  #get the start position of the gene in the chromosome
        chr_to= m[:end_chr].to_i #get the end podition of the gene in the chromosome
        
        entry.features do |feature| #consult all the features
          
          if feature.feature == 'exon' #and get the info of the exons only
            if feature.position =~ /^(?<from>\d+)..(?<to>\d+)/  #if the exon is in the forward strand (if it starts with the coordinates)
              coord= feature.position.match(/(?<from>\d+)..(?<to>\d+)/) #store it in the variable 'coord'
              
              unless exoncoordlist.include? coord #if it is not present, add to the list
                exoncoordlist << coord
              
                from= coord[:from].to_i #take the first and last coordinates of the exon
                to = coord[:to].to_i
                exon= seqembl.subseq(from, to) #and get the sequence between those coordinates
                exonlist << exon
                
                if exon.match(/(?=(cttctt))/i) #if it matchs with the motif 'cttctt'
                  geneinlist << "#{line}" #add the geneid to the list of genes that has the motif
                  
                  start0= exon.enum_for(:scan, /(?=(cttctt))/i).map { Regexp.last_match.begin(0)} #get the start position of the motif
                  start0.each {|e|#if there are more than one repetition, split them and consult each one
                    start= e + from #get the starting position in the gene
                    last= e + from + 5 #and the last
                    
                    #and same on the chromosome
                    chr_start= chr_from + start -1 
                    chr_end = chr_start + 5
                    
                    chr_coord = "#{chr_start}..#{chr_end}" #get the coordinates in chromosome
                    unless chrcoordlist.include? chr_coord
                      chrcoordlist << chr_coord
                    end
                    
                    gene_coordinates="#{start}..#{last}" # Get coordinates in gene
                    unless genecoordlist.include? gene_coordinates
                      genecoordlist << gene_coordinates                      
                    end
                
                    qual = feature.assoc 
                    attributelist << qual['note'] #store also the additional information of the exon
                    
                    if feature.locations.first.strand.to_s.match(/^1/) #it is always going to be in the forward if it is in this part of the
                      #conditional, but just in case, and to be easier the comprobation, lets do the whole 
                      strand= '+'
                      strandlist << strand
                    elsif feature.locations.first.strand.to_s.match(/^\-1/)
                      strand= '-'
                      strandlist << strand
                    end}
            
                end
              end
              
            elsif feature.position =~ /^complement\((?<from>\d+)..(?<to>\d+)\)/ #do the same with the reverse strand (the one that has complement
              #word before the coordinates)
              coord= feature.position.match(/^complement\((?<from>\d+)..(?<to>\d+)\)/)
              unless exoncoordlist.include? coord
                exoncoordlist << coord
              
                
                from= coord[:from].to_i
                to = coord[:to].to_i
                exon= seqembl.subseq(from, to)
                exonlist << exon
                
                if exon.match(/(?=(aagaag))/i)
                  geneinlist << "#{line}"
                     
                  start0= exon.enum_for(:scan, /(?=(aagaag))/i).map { Regexp.last_match.begin(0)}
                  start0.each {|e|
                    start= e + from
                    last= e + from + 5
                    chr_start= chr_from + start - 1
                    chr_end = chr_start + 5
                    
                    chr_coord = "#{chr_start}..#{chr_end}"
                    unless chrcoordlist.include? chr_coord
                      chrcoordlist << chr_coord
                    end                  
                    gene_coordinates="#{start}..#{last}" 
                    unless genecoordlist.include? gene_coordinates
                      genecoordlist << gene_coordinates
                      
                    end
                  
                    qual = feature.assoc
                    attributelist << qual['note']
                  
                    if feature.locations.first.strand.to_s.match(/^1/)
                      strand= '+'
                      strandlist << strand
                    elsif feature.locations.first.strand.to_s.match(/^\-1/)
                      strand= '-'
                      strandlist << strand
                    end}
                end
              end
            end          
          end

        end
      create_features(genecoordlist,emblseq,strandlist,'generepeats','cttctt', attributelist) #do the function to add the feature to the object
      #with the info in genes
      create_features(chrcoordlist,emblseq,strandlist,'chrmrepeats','cttctt', attributelist) #and in chromosomes
      
      emblseq.features.each do |feature| #check the features and stop when found 'generepeats'
        featuretype = feature.feature
        case
        when featuretype == "generepeats" 
          write_gff("gff_gene.gff",line,feature)#write a gff file with the gene info
        end
      end
      
      emblseq.features.each do |feature|#check the features and stop when found 'generepeats'
        featuretype = feature.feature
        case
        when featuretype == "chrmrepeats" 
          write_gff("gff_chr.gff",chr,feature) #write a gff file with the chromosome info
        end
      end         
       
      next if geneinlist.include? line #and if the gene isn't in the list of the genes that has the motif
      unless nocttgenes.include? line #unless it is already in the list of genes that doesn't have the motif
        nocttgenes << line #add to it to avoid repetitions
        report1= File.open("genes_without_cttctt.txt", "a+")
        report1.write("#{line}\n") #add the gene id
        report1.close
      end
      end  
    end
listfile.close