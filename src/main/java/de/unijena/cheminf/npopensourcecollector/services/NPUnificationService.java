package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProductRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
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

    private ArrayList<String> sourceNames;



    public void doWork(){

        System.out.println("NP unification InChi-key based");

        sourceNames = this.fetchSourceNames();

        System.out.println("SOURCES  "+sourceNames);

        List<String> uniqueInchiKeys = sourceNaturalProductRepository.findUniqueInchiKeys();

        for(String oinchikey: uniqueInchiKeys){


            String inchikey =  oinchikey.split("\"")[3];  //oinchikey.toString();
            System.out.println(inchikey);
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
            unp.found_in_databases = new HashSet<>();
            unp.absolute_smiles = new Hashtable<>();

            //associate the UniqueNaturalProduct entry to each of the sources
            for(SourceNaturalProduct snp : snpList){
                snp.setUniqueNaturalProduct(unp);
                sourceNaturalProductRepository.save(snp);

                //add annotations from SourceNaturalProducts

                //name
                //checking if name doesn't contain DB name
               boolean nameIsReal = true;
                for(String dbname: this.sourceNames){
                    if(snp.getName() != null && snp.getName().toLowerCase().contains(dbname)){
                        nameIsReal=false;
                    }
                }


                if(nameIsReal && snp.getName() != null && ((unp.getName() == null || unp.getName() =="") &&  snp.getName().length()>3)){


                        String name = snp.getName().trim();

                        String[] names = name.split("\\\n");


                        unp.setName(names[0]);
                        if (names.length > 1) {
                            for (int i = 1; i < names.length; i++) {
                                unp.synonyms.add(names[i]);
                            }
                        }


                }
                else if( unp.getName() != null && snp.getName() != null){
                    if(snp.getSource().equals("chebi")){

                        //replace name by ChebiName
                        unp.synonyms.add(unp.name);
                        unp.name = snp.getName().trim();
                    }
                    else {
                        unp.synonyms.add(snp.getName().trim());
                    }
                }

                //synonyms
                if(snp.getSynonyms() != null){

                    String[] synonyms;

                    for(String sy : snp.getSynonyms()){

                        String[] names = sy.split("\\\n");
                        for (int i = 0; i < names.length; i++) {
                            unp.synonyms.add(names[i].trim());
                        }
                    }

                }


                //species
                if(snp.organismText != null ){

                    unp.textTaxa.addAll(snp.organismText);
                }
                if(snp.taxid != null){
                    unp.taxid.addAll(snp.taxid);
                }
                if(unp.textTaxa.size()>1 && unp.textTaxa.contains("notax")){
                    unp.textTaxa.remove("notax");
                }


                //geo
                if(snp.getGeographicLocation() != null){
                    unp.geoLocation.addAll(snp.getGeographicLocation());
                }
                if(snp.getContinent() != null){
                    unp.geoLocation.add(snp.getContinent());
                }
                if(unp.geoLocation.size()>1 && unp.geoLocation.contains("nogeo")){
                    unp.geoLocation.remove("nogeo");
                }


                //refs
                if(snp.getCitation() != null){
                    for(String cit : snp.getCitation()){
                        unp.citationDOI.add(cit.trim());
                    }

                }

                //cas
                if(snp.getCas() != null && snp.getCas() != ""){
                    unp.setCas(snp.getCas() );
                }

                //database
                if(snp.getSource() != null){
                    unp.found_in_databases.add(snp.getSource());
                }

                //Absolute smiles (with stereochemistry)
                if(snp.getAbsoluteSmiles() != null && !snp.getAbsoluteSmiles().equals("")) {
                    if (unp.absolute_smiles.containsKey(snp.getAbsoluteSmiles())) {
                        unp.absolute_smiles.get(snp.getAbsoluteSmiles()).add(snp.getSource());

                    } else {
                        HashSet newSourceList = new HashSet();
                        newSourceList.add(snp.getSource());
                        unp.absolute_smiles.put(snp.getAbsoluteSmiles(), newSourceList);
                    }
                }
                else{
                    if(unp.absolute_smiles.containsKey("nostereo")) {
                        unp.absolute_smiles.get("nostereo").add(snp.getSource());
                    }
                    else{
                        HashSet newSourceList = new HashSet();
                        newSourceList.add(snp.getSource());
                        unp.absolute_smiles.put("nostereo", newSourceList);
                    }
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

        //TODO min cycle base (min number of rings)

        // count rings
        try {
            IRingSet rs = arf.findAllRings(im, 20);

            m.setMax_number_of_rings(rs.getAtomContainerCount());

            Cycles   cycles = Cycles.sssr(im);
            IRingSet rings  = cycles.toRingSet();
            m.setMin_number_of_rings(rings.getAtomContainerCount()); //SSSR


        } catch (CDKException e) {
            System.out.println("Too complex: "+m.getSmiles());
        }

        //compute molecular formula
        m.setMolecular_formula(mfm.getString(mfm.getMolecularFormula(im) ));

        //compute number of carbons, of nitrogens, of oxygens
        m.setNumber_of_carbons(mfm.getElementCount(mfm.getMolecularFormula(im), "C"));

        m.setNumber_of_oxygens(mfm.getElementCount(mfm.getMolecularFormula(im), "O"));

        m.setNumber_of_nitrogens(mfm.getElementCount(mfm.getMolecularFormula(im), "N"));

        m.setMolecular_weight( acm.getMass(im) );


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




    public ArrayList fetchSourceNames(){
        List<Object> uniqueSourceNames = sourceNaturalProductRepository.findUniqueSourceNames();
        ArrayList<String> tmpArray = new ArrayList<>();

        for(Object usn: uniqueSourceNames) {
            tmpArray.add(usn.toString().toLowerCase());
        }
        return tmpArray;
    }


}
