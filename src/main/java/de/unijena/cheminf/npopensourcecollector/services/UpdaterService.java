package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProductRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.List;
import java.util.Optional;

@Service
public class UpdaterService {

    @Autowired
    SourceNaturalProductRepository sourceNaturalProductRepository;

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    public void updateSourceNaturalProducts(){

        //get all sourceNaturalProduct

        List<SourceNaturalProduct> allSourceNaturalProducts = sourceNaturalProductRepository.findAll();

        for(SourceNaturalProduct snp : allSourceNaturalProducts){
            String unpid = snp.getUniqueNaturalProduct().getId();
            Optional<UniqueNaturalProduct> unp = uniqueNaturalProductRepository.findById(unpid);
            if(unp.isPresent()){
                UniqueNaturalProduct np = unp.get();
                snp.setUniqueNaturalProduct(np);
                sourceNaturalProductRepository.save(snp);
            }

        }
    }
}
