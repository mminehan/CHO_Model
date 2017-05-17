using LightXML
using JSON
using Requests


known_species_array = ParseXML()

print(known_species_array)












function ParseXML()
main = parse_file("cho_cobra.xml")

trunk = root(main)
println(name(trunk))
url_id = String[]
species_id = String[]
species_name = String[]

for c1 in child_nodes(trunk)  # c is an instance of XMLNode
    println(nodetype(c1))
    if is_elementnode(c1)
        branch = XMLElement(c1)  # this makes an XMLElement instance
        #println(name(branch))
        for c2 in child_nodes(branch)
          if is_elementnode(c2)
              branch2 = XMLElement(c2)  # this makes an XMLElement instance
              #println(name(branch2))
              if name(branch2) == "listOfSpecies"
                #println(nodetype(branch2))
                for c3 in child_nodes(branch2) #ITERATE THRU SPECIES HERE
                  if is_elementnode(c3)
                  species = XMLElement(c3)
                  #println(name(species))
                  if name(species) == "species"
                    for c4 in child_nodes(species)
                      if is_elementnode(c4)
                        annotation = XMLElement(c4)
                        #println(name(annotation))
                        if name(annotation) == "annotation"
                          for c5 in child_nodes(annotation)
                            if is_elementnode(c5)
                              RDF = XMLElement(c5)
                              #println(name(RDF))
                              if name(RDF) == "RDF"
                                for c6 in child_nodes(RDF)
                                  if is_elementnode(c6)
                                    Description = XMLElement(c6)
                                    #println(name(RDF))
                                    if name(Description) == "Description"
                                      #println(Description["is"])
                                      for c7 in child_nodes(Description)
                                        if is_elementnode(c7)
                                          is = XMLElement(c7)
                                          #println(name(is))
                                          if name(is) == "is"
                                            for c8 in child_nodes(is)
                                              if is_elementnode(c8)
                                                Bag = XMLElement(c8)
                                                #println(name(Bag))
                                                #println(Bag["li"])
                                                for c9 in child_nodes(Bag)
                                                  if is_elementnode(c9)
                                                    li = XMLElement(c9)
                                                    println(attribute(li,"resource"))
                                                    push!(url_id,attribute(li,"resource"))
                                                    push!(species_id,attribute(species,"id"))
                                                    push!(species_name,attribute(species,"name"))
                                                  end
                                                end
                                              end
                                            end
                                          end
                                        end
                                      end
                                    end
                                  end
                                end
                              end
                            end
                          end
                        end
                      end
                      end
                    end
                  end
                end
              end
          end
        end
    end
end

out = [species_name;species_id;url_id;]
return out

end
