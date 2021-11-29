require './geneinfo.rb'
require 'json'
require 'rest-client'
require './networks.rb'


# ruby maincode.rb prueba.txt 
Geneinfo.load_from_file(ARGV[0])

if File.file?("./networkreport.txt") 
  File.delete("./networkreport.txt") #If the file already exists, we delete it
end
if File.file?("./networkreport2.txt")
  File.delete("./networkreport2.txt") #If the file already exists, we delete it
end

Geneinfo.class_variable_get(:@@my_genes).each do |x| Geneinfo.annotate(x)
end
Geneinfo.class_variable_get(:@@my_genes).each do |y| Geneinfo.annotate_kegg(y)
end
Geneinfo.class_variable_get(:@@my_genes).each do |z| Geneinfo.annotate_interactors(z)
end
#Geneinfo.what_is_here
Networks.annotatethings
#Networks.what_is_here
Networks.class_variable_get(:@@allnetworks).each do |w| Networks.createNetworks(w)
end
Networks.class_variable_get(:@@allnetworks).each do |a| Networks.listinteractions(a)
end

puts "Conclusion: I think that I would need deeper networks to conclude something, but I think that there are no enough interactions between genes of the list to state that there is corregulation, so the paper is not very reliable"