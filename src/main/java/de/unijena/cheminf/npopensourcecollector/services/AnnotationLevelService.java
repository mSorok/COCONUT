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
        System.out.println("Evaluating annotations");
        List<String> allCoconutIds = uniqueNaturalProductRepository.findAllCoconutIds();

        for(String coconut_id : allCoconutIds){
            this.doWorkForOne(coconut_id);
        }
        System.out.println("done");
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

            np.textTaxa.remove("");
            if (!np.getTextTaxa().isEmpty() && !np.getTextTaxa().equals("notax") ){
                hasOrganism = true;
            }

            np.citationDOI.remove("");
            if (!np.citationDOI.isEmpty()) {
                hasLiterature = true;
            }

            if (!np.getFound_in_databases().isEmpty() && (np.found_in_databases.contains("chebi_np") || np.found_in_databases.contains("knapsack") || np.found_in_databases.contains("cmaup") || np.found_in_databases.contains("chembl_np") || np.found_in_databases.contains("npatlas"))) {
                hasTrustedSource = true;
            }

            boolean[] levels = {hasName, hasOrganism, hasLiterature, hasTrustedSource};

            annotationlevel = Booleans.countTrue(levels) + 1;

            np.annotationLevel = annotationlevel;

            uniqueNaturalProductRepository.save(np);
        }



    }
}
