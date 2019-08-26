package de.unijena.cheminf.npopensourcecollector.services;

import com.google.common.collect.Sets;
import de.unijena.cheminf.npopensourcecollector.mongocollections.NPSimilarity;
import de.unijena.cheminf.npopensourcecollector.mongocollections.NPSimilarityRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.similarity.DistanceMoment;
import org.openscience.cdk.similarity.Tanimoto;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.*;

@Service
public class SimilarityComputationService {

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    @Autowired
    NPSimilarityRepository npSimilarityRepository;

    @Autowired
    AtomContainerToUniqueNaturalProductService atomContainerToUniqueNaturalProductService;

    private Set<Set<UniqueNaturalProduct>> npPairs;

    public void computeSimilarities(){
        Fingerprinter fingerprinter = new Fingerprinter();

        System.out.println("Computing similarities");


        for(Set spair : npPairs){
            List<UniqueNaturalProduct> pair = new ArrayList<UniqueNaturalProduct>(spair);

            IAtomContainer mol1 = atomContainerToUniqueNaturalProductService.createAtomContainer(pair.get(0));
            IAtomContainer mol2 = atomContainerToUniqueNaturalProductService.createAtomContainer(pair.get(1));

            IBitFingerprint fingerprint1 = null;
            IBitFingerprint fingerprint2 = null;
            try {
                fingerprint1 = fingerprinter.getBitFingerprint(mol1);
                fingerprint2 = fingerprinter.getBitFingerprint(mol2);
                double tanimoto_coefficient = Tanimoto.calculate(fingerprint1, fingerprint2);

                if (tanimoto_coefficient>=0.5){

                    NPSimilarity newSimilarity = new NPSimilarity();
                    newSimilarity.setUniqueNaturalProduct1(pair.get(0));
                    newSimilarity.setUniqueNaturalProduct2(pair.get(1));
                    newSimilarity.setTanimoto(tanimoto_coefficient);
                    //newSimilarity.setDistanceMoment(distance_moment);

                    npSimilarityRepository.save(newSimilarity);
                }

            } catch (CDKException e) {
                e.printStackTrace();
            }

        }
        System.out.println("done");

    }


    public void generateAllPairs(){

        List<UniqueNaturalProduct> allNP = uniqueNaturalProductRepository.findAll();

        Set<UniqueNaturalProduct> npset = new HashSet<>(allNP);


        System.out.println("Computing pairs of NPs");

        this.npPairs = Sets.combinations( npset, 2);

        System.out.println("done");

    }
}
