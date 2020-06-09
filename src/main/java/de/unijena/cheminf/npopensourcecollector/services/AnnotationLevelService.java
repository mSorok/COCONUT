package de.unijena.cheminf.npopensourcecollector.services;

import com.google.common.primitives.Booleans;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.List;



@Service
public class AnnotationLevelService {

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    public void doWorkForAll(){
        List<String> allCoconutIds = uniqueNaturalProductRepository.findAllCoconutIds();

        for(String coconut_id : allCoconutIds){
            this.doWorkForOne(coconut_id);
        }

    }


    public void doWorkForOne(String coconut_id){

        UniqueNaturalProduct np = uniqueNaturalProductRepository.findByCoconut_id(coconut_id).get(0);

        if(np != null) {

            Integer annotationlevel = 1; //0 means it was never treated, 1 means that it was treated, but none of the

            boolean hasName = false;
            boolean hasOrganism = false;
            boolean hasLiterature = false;
            boolean hasTrustedSource = false;


            if (np.getName() != "" && np.getName() != null) {
                hasName = true;
            }

            if (!np.getTextTaxa().isEmpty()) {
                hasOrganism = true;
            }

            if (!np.citationDOI.isEmpty()) {
                hasLiterature = true;
            }

            if (!np.getFound_in_databases().isEmpty() && (np.found_in_databases.contains("chebi") || np.found_in_databases.contains("pubchem") || np.found_in_databases.contains("cmaup") || np.found_in_databases.contains("chembl"))) {
                hasTrustedSource = true;
            }

            boolean[] levels = {hasName, hasOrganism, hasLiterature, hasTrustedSource};

            annotationlevel = Booleans.countTrue(levels) + 1;

            np.annotationLevel = annotationlevel;

            uniqueNaturalProductRepository.save(np);
        }


    }
}
