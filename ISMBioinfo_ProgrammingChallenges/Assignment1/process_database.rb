require './seed_stock_data.rb'
require './gene_information.rb'
require './cross_data.rb'
require 'csv'


Gene.load_from_file(ARGV[0])
Seed.load_from_file(ARGV[1])
Cross_data.load_from_file(ARGV[2])

#puts Gene.how_many #some stuff to check if everything is ok
#puts Seed.how_many
#puts Cross_data.how_many
#puts Cross_data.what_is_here
#puts; puts
Seed.plant_grams
Seed.write_database(ARGV[3])
Cross_data.chisquare

