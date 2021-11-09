require './seed_stock_data.rb'
require './gene_information.rb'

class Cross_data

  @@number_of_crosses = 0 #start the set value, it is to know how many objects are created
  @@cross_obs = [] #start an array where are going to be all the objects
  attr_accessor :parent1
  attr_accessor :parent2
  attr_accessor :f2_Wild
  attr_accessor :f2_P1
  attr_accessor :f2_P2
  attr_accessor :f2_P1P2
  attr_accessor :linkedgene1 #this will be used in the def chisquare, to store the value of the gene name and make the final report
  attr_accessor :linkedgene2
    
  def initialize (params = {})
    @@number_of_crosses+= 1
    @parent1 = params.fetch(:parent1, '00000')
    @parent2 = params.fetch(:parent2, '00000')
    @f2_Wild = params.fetch(:f2_Wild, "000").to_f #convert to float because if not, the chisquare test gives a bad result
    @f2_P1 = params.fetch(:f2_P1, "00").to_f
    @f2_P2 = params.fetch(:f2_P2, "00").to_f
    @f2_P1P2 = params.fetch(:f2_P1P2, "00").to_f
    @linkedgene1 = params.fetch(:linkedgene1, "Unknowngene1")
    @linkedgene2 = params.fetch(:linkedgene2, "Unknowngene2")
    @@cross_obs << self # 'self' it’s a Ruby keyword that gives you access to the current object, so I add each object to the array
  end
    
  def Cross_data.load_from_file(file_name)
    myfile = File.open(file_name) #I tried to do File.open(file_name, headers: true) but didn't work, i found that .drop(n) in readlines skips n lines from the begining
    File.readlines(myfile).drop(1).each do |line|
      p1, p2, f2_wild, f2_p1, f2_p2, f2_p1p2 = line.split("\t") #stores in each variable the value of each column separated with tabs
      self.new(parent1: p1, parent2: p2, f2_Wild: f2_wild, f2_P1: f2_p1, f2_P2: f2_p2, f2_P1P2: f2_p1p2) #create new objects assigning the values of the previous variables to the attributes 
    myfile.close
    end
  end
    
  def Cross_data.how_many
    return @@number_of_crosses #how many objects are created, just for my information while doing the code
  end
    
  def Cross_data.what_is_here
    return @@cross_obs.each {|x| puts x.parent1}#to know if the objects are being well created while doing the code (i have had some struggles with this)
  end
    
  def Cross_data.chisquare  #voy a probar con sum() sino hay que escribirla entera
    @@cross_obs.each {|x|
      all= (x.f2_Wild + x.f2_P1 + x.f2_P2 + x.f2_P1P2)
      test = ((x.f2_Wild - (all*9/16))**2/(all*9/16)) + ((x.f2_P1 - (all*3/16))**2/(all*3/16)) + ((x.f2_P2 - (all*3/16))**2/(all*3/16)) + ((x.f2_P1P2 - (all*1/16))**2/(all*1/16))
      if test > 7.81 #as we have 4 phenotypes, the degrees of freedom will be 3 and I took a 0.05 alfa value
        Seed.class_variable_get(:@@stock_obs).each {|z| #for the cases that H1: being genetically linked is acepted, consult the seed class obkects and,
          if x.parent1 == z.seed_stock #if the mutation of the parent1 matches the mutation of the seed_stock of Seed.class,
            Gene.class_variable_get(:@@gene_obs).each {|y|  #consult the Gene class objects and,
              if y.gene_ID == z.mutant_Gene_ID #then, if the mutant_gene_id of Seed class match the gene_id of the Gene class,
                  x.linkedgene1 = y.gene_name #store in the "linkedgene1" cross' attribute the value of the gene_name that is linked, so that it is possible to acces that value in the final report "genename is genetically..."
                  y.linked_to = x.linkedgene1 #add the information to which gene is genetically linked as a value of the attribute @linked_to
              end}
          elsif x.parent2 == z.seed_stock #do the same here
            Gene.class_variable_get(:@@gene_obs).each {|y|
              if y.gene_ID == z.mutant_Gene_ID
                  x.linkedgene2 = y.gene_name
                  y.linked_to = x.linkedgene2
              end}
          puts "#{x.linkedgene1} is genetically linked to #{x.linkedgene2} with chisquare score #{test}"
          end}   
      end}
  end
end