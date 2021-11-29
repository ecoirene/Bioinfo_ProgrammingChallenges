
class Seed

  @@number_of_seeds = 0 #start the set value, it is to know how many objects are created
  @@stock_obs = [] #start an array where are going to be all the objects
  attr_accessor :seed_stock
  attr_accessor :mutant_Gene_ID
  attr_accessor :last_Planted
  attr_accessor :storage
  attr_accessor :grams_Remaining
  
  def initialize (params = {})
    @@number_of_seeds+=1  # increment the number of seeds by one
    @seed_stock = params.fetch(:seed_stock, "00000")
    @mutant_Gene_ID = params.fetch(:mutant_Gene_ID, 'AT0G00000')
    @last_Planted = params.fetch(:last_Planted, '00/00/0000')
    @storage = params.fetch(:storage, 'cama00')
    @grams_Remaining = params.fetch(:grams_Remaining, "00").to_i
    @@stock_obs << self # 'self' it’s a Ruby keyword that gives you access to the current object, so I add each object to the array
  end
  
  def Seed.load_from_file(file_name)
    myfile = File.open(file_name)#I tried to do File.open(file_name, headers: true) but didn't work, i found that .drop(n) in readlines skips n lines from the begining
    File.readlines(myfile).drop(1).each do |line|
      seed_stock, mutant_gene_ID, last_planted, storage, grams_remaining = line.split("\t")  #stores in each variable the value of each column separated with tabs
      self.new(seed_stock: seed_stock, mutant_Gene_ID: mutant_gene_ID, last_Planted: last_planted, storage: storage, grams_Remaining: grams_remaining) #create new objects assigning the values of the previous variables to the attributes 
    myfile.close
    end
  end

  def Seed.how_many
    return @@number_of_seeds #how many objects are created, just for my information while doing the code
  end
  
  def Seed.what_is_here
    return @@stock_obs.each {|x| puts x.grams_Remaining} #to know if the objects are being well created while doing the code (i have had some struggles with this)
  end
  
  def Seed.plant_grams
    @@stock_obs.each {|x|
      if x.grams_Remaining == 7
        puts "WARNING: we have just run out of seeds in #{x.seed_stock}" #if the value of that attribute is equal to 7, will take those but print a warning message
        x.grams_Remaining = x.grams_Remaining - 7
        x.last_Planted = Time.new.strftime("%d/%m/%Y") #and update the date
      elsif x.grams_Remaining > 7 #if there is more than 7, we can just take them without problem
        x.grams_Remaining = x.grams_Remaining - 7
        x.last_Planted = Time.new.strftime("%d/%m/%Y")
      else
        puts "WARNING: there are no enough seeds to plant in #{x.seed_stock}" #if it is less than 7, will no take anything and will print a warning message
      end}
  end
  
  def Seed.write_database(new_file)
    updatedfile=File.open(new_file, 'w') #create a new database
    updatedfile << "Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining\n" #add the titles of the columns
    @@stock_obs.each {|row|
      updatedfile.write("#{row.seed_stock}\t#{row.mutant_Gene_ID}\t#{row.last_Planted}\t#{row.storage}\t#{row.grams_Remaining}\n")# add the values of each object to each row
      }
    end
  end