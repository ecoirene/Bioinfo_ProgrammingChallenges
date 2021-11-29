
class Gene

    @@number_of_genes = 0 #start the set value, it is to know how many objects are created
    @@gene_obs = [] #start an array where are going to be all the objects
    attr_accessor :gene_ID
    attr_accessor :gene_name
    attr_accessor :mutant_phenotype
    attr_accessor :linked_to
    
    def initialize (params = {})
        @@number_of_genes+=1 
        @gene_ID = params.fetch(:gene_ID, 'AT0G00000')
        @gene_name = params.fetch(:gene_name, 'SomeGene')
        @mutant_phenotype = params.fetch(:mutant_phenotype, "Some details")
        @linked_to = params.fetch(:linked_to, "Unknown gene")
        @@gene_obs << self # 'self' it’s a Ruby keyword that gives you access to the current object, so I add each object to the array
    end
    
    def Gene.load_from_file(file_name) 
        myfile = File.open(file_name) #I tried to do File.open(file_name, headers: true) but didn't work, i found that .drop(n) in readlines skips n lines from the begining
        File.readlines(myfile).drop(1).each do |line|
          gene_ID, gene_name, mutant_phenotype = line.split("\t") #stores in each variable the value of each column separated with tabs
          id_checking = Regexp.new(/A[Tt]\d[Gg]\d\d\d\d\d/) #check the id of arabidopsis
          unless id_checking.match(gene_ID)
            abort("WARNING, gene identifier #{id_checking} has not the right format")
          end
          self.new(gene_ID: gene_ID, gene_name: gene_name, mutant_phenotype: mutant_phenotype) #create new objects assigning the values of the previous variables to the attributes 
        myfile.close
        end
    end

    def Gene.how_many
        return @@number_of_genes #how many objects are created, just for my information while doing the code
    end
    
end