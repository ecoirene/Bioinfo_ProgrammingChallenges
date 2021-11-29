require './geneinfo.rb'
require 'rest-client'
require 'json'

#sorry for the mess, I will try to be cleaner and more efficient in the future

def fetch(url, headers = {accept: "*/*"}, user = "", pass="") 
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

class Networks
  
  @@allnetworks = []
  attr_accessor :begining #first gene of the network, the same as @gene_id in the Geneinfo class
  attr_accessor :intact2 #  interactors of the input list (is the same as @intact in Geneinfo class)
  attr_accessor :network_prots #Array with all the genes that interacts with intact2
  
  def initialize (params = {})
    @begining = params.fetch(:begining, 'AT0G00000')
    @network_prots = params.fetch(:network_prots, [])
    @intact2 = params.fetch(:intact2, []) 
    @@allnetworks << self
    
  end

  def Networks.createNetworks(w) #search the interactions of the interactors
      address = fetch("http://bar.utoronto.ca:9090/psicquic/webservices/current/search/interactor/#{w.intact2}?format=tab25")
      
      if address
        data = address.to_s.split("\n") #this block is the same as the one used in geneinfo class
        data.each do |line|  
        element = line.split("\t")   
        gene1 = element[2].sub(/tair:/, "").upcase  
        gene2 = element[3].sub(/tair:/, "").upcase  
        taxid1 = element[9].sub(/taxid:/, "").to_i
        taxid2 = element[10].sub(/taxid:/, "").to_i
        score = element[14].sub(/intact-miscore:/, "").to_f
        #puts gene1, gene2, score
        if taxid1 != 3702 or taxid2 != 3702
          next
        elsif score < 0.45
          next
        elsif gene1 == w.begining or gene2 == w.begining #if one of them is the same as the first gene of the network (just the inverse interaction), next
          next
        else #in the case that they are new interactions and have good quality
          case
          when gene1 == w.intact2 #when gene 1 is the query id, store the gene 2 in the array
            w.network_prots << gene2 
          else #and viceversa
            w.network_prots << gene1
          end
        end
        end  
      end
      
    w.network_prots.each {|f| #open the array of the last interactors for each object
      Geneinfo.class_variable_get(:@@my_genes).each {|g| #and have access to the objects of the other class
      if g.gene_id == f #if the input gene_id and the last interactor is the same
        #puts "#{w.begining} interacts with #{w.intact2} which also interacts with #{f}"
        myfile = File.open("./networkreport.txt", "a+")
        #write the gene_id (synonimus of w.begining) > the 2nd interactor* > and finally the last, which interacts with the 2nd and is part of the input list
        #*I have created as many objects as 2nd interactors are, so w.intact2 is finally not an array
        myfile.puts "#{w.begining} interacts with #{w.intact2} which also interacts with #{f}\n"
        myfile.close
      end}}
  end


  def Networks.annotatethings #import the needed info to this class objects
    Geneinfo.class_variable_get(:@@my_genes).each {|x|
      x.intact.each {|inter|
        Networks.new(intact2: inter, begining: x.gene_id)}}
  end
  
  def Networks.what_is_here #just for me to try to understand what I was doing
    return Networks.class_variable_get(:@@allnetworks).each {|x| puts "#{x.begining} interacts with #{x.intact2}"}
  end
  
  def Networks.listinteractions(a) #create another report with the genes of the input list that interact with each other
      Geneinfo.class_variable_get(:@@my_genes).each {|c|
        if c.gene_id == a.intact2  #if the interactor is in the list
          unless a.begining == a.intact2 #but is not the same of the query gene, write that info in the txt file
            #puts "The genes of the list #{a.begining} and #{a.intact2} interact with each other"
            myfile1 = File.open("./networkreport2.txt", "a+")
            myfile1.puts "The genes of the list #{a.begining} and #{a.intact2} interact with each other \n"       
            myfile1.close
          end
        end}
  end
  
end

