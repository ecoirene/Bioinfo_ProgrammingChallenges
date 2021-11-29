require 'rest-client'
require 'json'


class Geneinfo

  @@my_genes = []              #Array with all the genes
  attr_accessor :gene_id       #Gene id
  attr_accessor :accession     #Array of accessions
  attr_accessor :kegg          #Array with kegg reaction ids
  attr_accessor :kegg_pathway  #Hash kegg reaction id => pathway
  attr_accessor :go            #Hash go id => go term
  attr_accessor :intact        #Array with the interactors of each gene
  
  
  
  def initialize (params = {})
    
    @gene_id = params.fetch(:gene_id, 'Unknown_id')
    @accession = []
    @go = {}
    @kegg = []
    @intact = Array.new
    @kegg_pathway = {}
    
    @@my_genes << self #store the objects in the array
    
  end
  
  def Geneinfo.load_from_file(filename) #enter the agis of the input list in the gene_id attribute 
    myfile = File.open(filename)
    File.readlines(myfile).each {|line|  
      self.new(gene_id:line.upcase.strip)} 
      myfile.close
    end
  end
  
  def fetch(url, headers = {accept: "*/*"}, user = "", pass="") #to handle errors
    response = RestClient::Request.execute({
      method: :get,
      url: url.to_s,
      user: user,
      password: pass,
      headers: headers})
    return response
    rescue RestClient::ExceptionWithResponse => e #si pasa este error, rescatarlo con ese codigo (cada error con su rescue)
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
  
  def Geneinfo.annotate(x)
    #Using the gene id annotates other info (accessions, KEGG, GO terms)
      puts "http://togows.dbcls.jp/entry/uniprot/#{x.gene_id}.json"
      address = fetch("http://togows.dbcls.jp/entry/uniprot/#{x.gene_id}.json")
      
      if address #if there are no errors
        data = JSON.parse(address.body) #take the info
        #puts data
        x.accession = data[0]["accessions"] #insert the info in the attribute of each object
        #puts x.accession
        x.kegg =  data[0]["dr"]["KEGG"] #insert the info in the attribute of each object      
        #puts x.kegg
        #GO
        unless data[0]["dr"]["GO"].nil?
          data[0]["dr"]["GO"].each do |annot| #insert the info in the attribute of each object, sometimes, there is no GO term
            #puts annot
            (id, ont)=annot
            #puts ont if ont.match(/^P:/)
            if ont.match(/^P:/) #only those that have biologic process, related to the P
              x.go[id] = ont  #Create a hash with go id => go term 
            end
          end
        end
      else
        puts "the Web call failed - see STDERR for details..."    #if there is any error, print this warning mesage  
        
      end

  end 
  def Geneinfo.annotate_interactors(z) # the utoronto db is much better than the togows for taking the interactors
    address = fetch("http://bar.utoronto.ca:9090/psicquic/webservices/current/search/interactor/#{z.gene_id}?format=tab25")
    
    if address
      data = address.to_s.split("\n") #each line has info of 1 interaction
      data.each do |line|  
        element = line.split("\t")   #separate the columns in elements
        gene1 = element[2].sub(/tair:/, "").upcase  #I only want the id so I delete the  'tair:'
        gene2 = element[3].sub(/tair:/, "").upcase  
        taxid1 = element[9].sub(/taxid:/, "").to_i #the taxid is useful to filter by species, since I know that the Arabidopsis taxid is 3702
        taxid2 = element[10].sub(/taxid:/, "").to_i
        score = element[14].sub(/intact-miscore:/, "").to_f #medium quality starts in 0.45
        #puts gene1, gene2, 
        if taxid1 != 3702 or taxid2 != 3702
          next
        elsif score < 0.45
          next
        else
          case
          when gene1 == z.gene_id #if the input gene is the first that appear, select the 2nd and store in the array @intact
            z.intact << gene2
            
          else 
            z.intact << gene1 #and viceversa           
          end
       
        end
      end
    end
      
  end
    
  def Geneinfo.annotate_kegg(y) #the kegg pathways are in other part of the db but we can acces them this way
    address = fetch("http://togows.dbcls.jp/entry/kegg-reaction/ath:#{y.gene_id}.json") 
    if address
      data = JSON.parse(address.body)
      y.kegg_pathway = data[0]['pathways']
    else
      puts "the Web call failed - see STDERR for details..."
    end
  end
  
  def Geneinfo.what_is_here #just to consult while coding
    return Geneinfo.class_variable_get(:@@my_genes).each {|x| puts x.intact}
  end
  
  