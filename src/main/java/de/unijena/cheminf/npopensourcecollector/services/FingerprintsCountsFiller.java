package de.unijena.cheminf.npopensourcecollector.services;


import de.unijena.cheminf.npopensourcecollector.mongocollections.PubFingerprintsCounts;
import de.unijena.cheminf.npopensourcecollector.mongocollections.PubFingerprintsCountsRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.Hashtable;
import java.util.List;

@Service
public class FingerprintsCountsFiller {

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    @Autowired
    PubFingerprintsCountsRepository pubFingerprintsCountsRepository;


    public void doWork(){

        System.out.println("Generating PubChemFingerprints counts");

        Hashtable<Integer, Integer> counts = new Hashtable<>();


        List<String> allCoconutIds = uniqueNaturalProductRepository.findAllCoconutIds();

        for(String coconut_id : allCoconutIds){
            UniqueNaturalProduct np = uniqueNaturalProductRepository.findByCoconut_id(coconut_id).get(0);

            if(!np.pubchemFingerprint.isEmpty()){
                for(Integer bitOnIndex : np.pubchemFingerprint){

                    if(counts.containsKey(bitOnIndex)){
                        //add count
                        counts.put(bitOnIndex, counts.get(bitOnIndex) + 1);

                    }else{
                        //create new key
                        counts.put(bitOnIndex, 1);
                    }

                }
            }

        }

        for(Integer bitOnIndex: counts.keySet() ){
            PubFingerprintsCounts pubFingerprintsCounts = new PubFingerprintsCounts();
            pubFingerprintsCounts.id = bitOnIndex;
            pubFingerprintsCounts.count = counts.get(bitOnIndex);
            pubFingerprintsCountsRepository.save(pubFingerprintsCounts);

        }




        System.out.println("done");

    }



}
