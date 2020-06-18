package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

@Service
public class NamingService {

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    public void namesToLowcase(){

        System.out.println("Converting names to low case");
        List<String> allCoconutIds = uniqueNaturalProductRepository.findAllCoconutIds();

        for(String coconut_id : allCoconutIds){
            this.npNameToLowcase(coconut_id);
        }
        System.out.println("done");

    }


    public void npNameToLowcase(String coconut_id){

        UniqueNaturalProduct np = uniqueNaturalProductRepository.findByCoconut_id(coconut_id).get(0);

        if(np != null) {

            if(np.name != null && np.name != ""){

                np.name = np.name.toLowerCase();

                if(!np.synonyms.isEmpty()){
                    List<String> lowSynonyms = np.synonyms.stream()
                            .map(String::toLowerCase)
                            .collect(Collectors.toList());
                    np.synonyms = new HashSet<>(lowSynonyms);
                }

            }
            uniqueNaturalProductRepository.save(np);
        }

    }
}
