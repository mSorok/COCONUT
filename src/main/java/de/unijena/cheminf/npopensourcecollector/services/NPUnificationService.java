package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProductRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.function.IntBinaryOperator;

@Service
public class NPUnificationService {

    @Autowired
    SourceNaturalProductRepository sourceNaturalProductRepository;

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    @Autowired
    AtomContainerToUniqueNaturalProductService atomContainerToUniqueNaturalProductService;

    public void doWork(){

        System.out.println("NP unification InChi-key based");

        List<Object> uniqueInchiKeys = sourceNaturalProductRepository.findUniqueInchiKeys();

        for(Object oinchikey: uniqueInchiKeys){


            String inchikey = oinchikey.toString();
            //System.out.println(inchikey);
            List<SourceNaturalProduct> snpList = sourceNaturalProductRepository.findBySimpleInchiKey(inchikey);
            //System.out.println(snpList.get(0).simpleInchiKey+" "+snpList.get(0).source);

            //create a new UniqueNatural product

            UniqueNaturalProduct unp = new UniqueNaturalProduct();

            unp.setInchikey(inchikey);
            unp.setInchi(snpList.get(0).simpleInchi);
            unp.setSmiles(snpList.get(0).simpleSmiles);
            unp.setTotal_atom_number(snpList.get(0).getTotalAtomNumber());
            unp.setHeavy_atom_number(snpList.get(0).getHeavyAtomNumber());

            unp = uniqueNaturalProductRepository.save(unp);

            unp.synonyms = new HashSet<>();
            unp.textTaxa = new HashSet<>();
            unp.taxid = new HashSet<>();
            unp.geoLocation = new HashSet<>();
            unp.citationDOI = new HashSet<>();

            //associate the UniqueNaturalProduct entry to each of the sources
            for(SourceNaturalProduct snp : snpList){
                snp.setUniqueNaturalProduct(unp);
                sourceNaturalProductRepository.save(snp);

                //add annotations from SourceNaturalProducts

                //name
                if((unp.getName() == null || unp.getName() =="") && snp.getName() != null){

                    String name = snp.getName().trim();

                    String [] names = name.split("\\\n");


                    unp.setName(names[0]);
                    if(names.length>1){
                        for(int i=1; i<names.length;i++){
                            unp.synonyms.add(names[i]);
                        }
                    }
                }
                else if( unp.getName() != null && snp.getName() != null){
                    unp.synonyms.add(snp.getName().trim());
                }

                //synonyms
                if(snp.getSynonyms() != null){
                    unp.synonyms.addAll(snp.getSynonyms());
                }


                //species
                if(snp.organismText != null ){
                    unp.textTaxa.addAll(snp.organismText);
                }
                if(snp.taxid != null){
                    unp.taxid.addAll(snp.taxid);
                }


                //geo
                if(snp.getGeographicLocation() != null){
                    unp.geoLocation.addAll(snp.getGeographicLocation());
                }
                if(snp.getContinent() != null){
                    unp.geoLocation.add(snp.getContinent());
                }


                //refs
                if(snp.getCitation() != null){
                    unp.citationDOI.addAll(snp.getCitation());
                }
            }

            unp = uniqueNaturalProductRepository.save(unp);

            //compute molecular parameters for the UniqueNaturalProduct

            unp = computeAdditionalMolecularFeatures(unp);
            uniqueNaturalProductRepository.save(unp);

        }



    }









    public UniqueNaturalProduct computeAdditionalMolecularFeatures(UniqueNaturalProduct m){
        AllRingsFinder arf = new AllRingsFinder();
        MolecularFormulaManipulator mfm = new MolecularFormulaManipulator();
        AtomContainerManipulator acm = new AtomContainerManipulator();

        IAtomContainer im = atomContainerToUniqueNaturalProductService.createAtomContainer(m);

        // count rings
        try {
            IRingSet rs = arf.findAllRings(im, 15);

            m.setNumber_of_rings(rs.getAtomContainerCount());


        } catch (CDKException e) {
            System.out.println("Too complex: "+m.getSmiles());
        }

        //compute molecular formula
        m.setMolecular_formula(mfm.getString(mfm.getMolecularFormula(im) ));

        //compute number of carbons, of nitrogens, of oxygens
        m.setNumber_of_carbons(mfm.getElementCount(mfm.getMolecularFormula(im), "C"));

        m.setNumber_of_oxygens(mfm.getElementCount(mfm.getMolecularFormula(im), "O"));

        m.setNumber_of_nitrogens(mfm.getElementCount(mfm.getMolecularFormula(im), "N"));

        m.setMolecular_weight( acm.getMolecularWeight(im) );


        // cleaning the NaNs
        if(m.getMolecular_weight().isNaN()){
            m.setMolecular_weight(0.0);
        }


        //get bond count
        IBond[] bonds = acm.getBondArray(im);
        int bondCount = 0;
        for(IBond b : bonds){
            if(b.getAtomCount() == 2){
                if(!b.getAtom(0).getSymbol().equals("H") && !b.getAtom(1).getSymbol().equals("H")){
                    bondCount++;
                }
            }
        }
        m.setBond_count(bondCount);
        return(m);
    }
}
